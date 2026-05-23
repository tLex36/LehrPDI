import requests
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import random
import xml.etree.ElementTree as ET
import json
import time
import zipfile
import io
import openpyxl
import pandas as pd
import re
import os
from scipy.stats import binom 
# import torch
# import torch.nn as nn
import matplotlib.pyplot as plt
import numpy as np
#from aiupred import AIUPred


NCBI_API_KEY = "1df5ddb36d910c631f330dea9da56ed22809"  # or "your_key_here"
ENSEMBL_VEP_PAYLOAD_SIZE = 200
aa_map = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    # Add more if needed
}
# class AIUPredMPS(AIUPred):
#    def _setup_device(self, force_cpu, gpu_num):
#        if torch.backends.mps.is_available() and not force_cpu:
#            self.device = torch.device("mps")
#        elif torch.cuda.is_available() and not force_cpu:
#            self.device = torch.device(f"cuda:{gpu_num}")
#        else:
#            self.device = torch.device("cpu")
def get_mane_transcript(gene_symbol):
    mane_df = pd.read_csv("mane_transcripts.csv")
    row = mane_df[mane_df["symbol"] == gene_symbol]
    return row.iloc[0]["Ensembl_nuc"].split(".")[0]#[["Ensembl_nuc"]].to_dict("records")
def get_ensemblVEP_aa_consequences(variants):

    payload = create_ensemblVEP_payload(variants)
    total_ensemblVEP_consequences = [get_amino_acid_consequences_batch(i) for i in payload] #v stores all the ensemblVEP predicted consequences for each variant called for. That includes the transcript consequence, predicted sift score, polyphen prediction, and transcript_ID
    total_ensemblVEP_consequences = [item for sublist in total_ensemblVEP_consequences for item in sublist]#flatten the batches into 1 long list
    return total_ensemblVEP_consequences
    #Now each element of the list is 1 variant, with it's ensemblVEP data packed into it. ie. v[0] = all the VEP info
def create_ensemblVEP_payload(variants):
    payload_size = ENSEMBL_VEP_PAYLOAD_SIZE
    #Pass in the list of grch38_pos variants
    #Returns the list of payloads to iterate through in 200 variant chunks

    #convert the variant strings from #-######-ref-alt format to #:#####-#####/alt
    variants = [v.split("-")[0]+":"+v.split("-")[1]+"-"+v.split("-")[1]+"/"+v.split("-")[3] for v in variants]
    #split the variants into as many 200/payload_size chunked lists so that I can pass these onto the get_amino_acid_consequences_batch() function
    payload = [variants[i:i+payload_size]for i in range(0,len(variants),payload_size)]
    return payload
def get_amino_acid_consequences_batch(variants): #pass in a list of up to 200 variants and get that list back
    url = "https://rest.ensembl.org/vep/human/region"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }


    payload = {"variants": variants, "LoF":1}  # list of "chrom:pos-pos/alt" strings
    response = requests.post(url, headers=headers, json=payload)
    response.raise_for_status()
    return response.json()
def my_parse_variants(variation):
    v = {}
    v["clinvar_id"]   = variation.get("VariationID", "N/A")
    v["name"]         = variation.get("VariationName", "N/A")
    v["accession"]    = variation.get("Accession", "N/A")
    #I need to find the protein consequence of the missense variant and store that as a protein_location, protein_ref, and protein_alt 
    


    interp = variation.find(".//Interpretation")
    if interp is not None:
        desc = interp.find("Description")
        v["clinical_significance"] = desc.text if desc is not None else "N/A"
    else:
        v["clinical_significance"] = "N/A"

    consequences = list({
        mc.get("Type", "")
        for mc in variation.findall(".//MolecularConsequence")
        if mc.get("Type")
    })
    v["molecular_consequences"] = consequences

    location = variation.find(
        ".//SimpleAllele/Location/SequenceLocation[@Assembly='GRCh38']"
    )
    if location is None:
        for loc in variation.findall(".//SequenceLocation[@Assembly='GRCh38']"):
            if loc.get("referenceAlleleVCF") is not None:
                location = loc
                break

    if location is not None:
        v["chromosome"]     = location.get("Chr", "N/A")
        v["position_start"] = location.get("start", "N/A")
        v["position_stop"]  = location.get("stop", "N/A")
        v["ref_allele"]     = location.get("referenceAlleleVCF", "N/A")
        v["alt_allele"]     = location.get("alternateAlleleVCF", "N/A")
        v['variant_id']     = "-".join([v['chromosome'],v['position_start'],v['ref_allele'],v['alt_allele']])
    else:
        v["chromosome"] = v["position_start"] = v["position_stop"] = "N/A"
        v["ref_allele"] = v["alt_allele"] = "N/A"

    #I want to extract from the name field the protein consequence of the missense variant and store that as a protein_location, protein_ref, and protein_alt
    #The name field is usually in the format "NM_000123.4:c.123A>G (p.Lys41Arg)" or similar. I can use a regex to extract the protein change information from the parentheses. The protein change is usually in the format "p.AAA123BBB" where AAA is the reference amino acid, 123 is the position, and BBB is the alternate amino acid. I can use a regex like r"\(p\.([A-Za-z]+)(\d+)([A-Za-z]+)\)" to extract these components.
    name = v["name"] or ""

    v["protein_ref"], v["protein_location"], v["protein_alt"] = parse_protein_consequence(name)
    review = variation.find(".//ReviewStatus")
    v["review_status"] = review.text if review is not None else "N/A"

    xref = variation.find(".//XRef[@DB='dbSNP']")
    v["rsid"] = f"rs{xref.get('ID')}" if xref is not None else "N/A"

    return v
def parse_protein_consequence(name):
    match = re.search(r"p\.([A-Za-z]+)(\d+)([A-Za-z]+)", name)
    if match:
        ref = match.group(1)
        loc = int(match.group(2))  # convert position to integer
        alt = match.group(3)
        #I want to convert the three letter amino acid code to the single letter amino acid code. I can use a dictionary to map the three letter codes to single letter codes. For example, {"Ala": "A", "Arg": "R", "Asn": "N", ...}
        aa_map = {
            "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
            "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
            "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
            "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
            # Add more if needed
        }
        ref = aa_map.get(ref, ref)
        alt = aa_map.get(alt, alt)
    else: ref = loc = alt = "N/A"
    return ref, loc, alt
def get_cleaned_truncation_variants_list(tmp_collected_trunc_var_VEP, transcript_id):
    tmp_trunc_vars = []
    results = []
    early_trunc_terms = ['frameshift_variant','protein_altering_variant','stop_gained']
    for j in range(len(tmp_collected_trunc_var_VEP)): 
            if 'transcript_consequences' in tmp_collected_trunc_var_VEP[j]:
                #I next need to check that my MANE transcript is a truncation
                #print(""+str([i for i in tmp_collected_trunc_var_VEP[j]['transcript_consequences'] if i['transcript_id'] == grin2b_mane][0]['amino_acids'])+" at index: "+str(j))

                print(f"set test: {[i for i in tmp_collected_trunc_var_VEP[j]['transcript_consequences'] if i['transcript_id'] == transcript_id]}")

                if (bool(set([i['consequence_terms'] for i in tmp_collected_trunc_var_VEP[j]['transcript_consequences'] if i['transcript_id'] == transcript_id][0]) & set(early_trunc_terms))):# checks if MANE consequence of variant is actually an early truncation
                    aa_tmp = str([i for i in tmp_collected_trunc_var_VEP[j]['transcript_consequences'] if i['transcript_id'] == transcript_id][0]['amino_acids'])
                    loc_tmp = str([i for i in tmp_collected_trunc_var_VEP[j]['transcript_consequences'] if i['transcript_id'] == transcript_id][0]['protein_start'])
                    if not aa_tmp[-1] in aa_map.values():
                        tmp_trunc_vars.append(tmp_collected_trunc_var_VEP[j])
                        print("truncating variant of aa/index is : "+aa_tmp+" "+str(loc_tmp)+"    truncation found at index: "+str(j))
    print(len(tmp_trunc_vars))
    return tmp_trunc_vars
#Get truncating variants from ClinVar and GnomAD and returns them as lists with transcript properties
def collate_truncation_variants(gene, transcript_id = None):
    if transcript_id is None:
        transcript_id = get_mane_transcript(gene)
    #Fetch all clinvar variants
    gene_root, gene_tmp_var = fetch_all_clinvar_for_gene(gene,api_key=NCBI_API_KEY)
    variations = gene_root.findall(".//VariationArchive")
    print(f'fetched variants from clinvar {len(variations)}')
    gene_all_clinvar = [my_parse_variants(v) for v in variations]
    #collected_trunc_var = grin2b_all_clinvar['frameshift_variants'] + grin2b_all_clinvar['nonsense_variants']
    collected_clinvar_ids = list(set(
        [i['variant_id']for i in gene_all_clinvar if 'variant_id' in i]
        )) #list(set()) is to get rid of any redundant variants that are listed as both frameshift and nonsense
    gene_trunc_clinvar_VEP = get_ensemblVEP_aa_consequences(collected_clinvar_ids)
    print(f"fetched transcription consequences from ensemblVEP for truncation variants from clinvar {len(gene_trunc_clinvar_VEP)}")
    gene_trunc_clinvar_list = get_cleaned_truncation_variants_list(gene_trunc_clinvar_VEP, transcript_id)

    #Fetch all gnomad variants
    gene_gnomad = fetch_all_gnomad_variants(gene)

    gnomad_early_trunc_terms = ['frameshift_variant','protein_altering_variant','stop_gained',]
    gene_early_trunc_gnomad = [i for i in gene_gnomad if i['consequence'] in gnomad_early_trunc_terms]
    collected_gnomad_ids = list(set([i['variant_id'] for i in gene_early_trunc_gnomad])) #Collects the list of variant_ids in gnomad so that it can be passed to ensemblVEP
    gene_trunc_gnomad_VEP = get_ensemblVEP_aa_consequences(collected_gnomad_ids)
    gene_trunc_gnomad_list = get_cleaned_truncation_variants_list(gene_trunc_gnomad_VEP, transcript_id)

    return gene_trunc_clinvar_list, gene_trunc_gnomad_list
def fetch_all_gnomad_variants(gene_symbol):
    #This function will query gnomAD for missense variants in the given gene symbol and return a list of missense variants with their positions and allele frequencies
    #I can use the gnomAD API to fetch this information. The endpoint for fetching variants by gene is "https://gnomad.broadinstitute.org/api/variants/search?query=gene:{gene_symbol}&variantType=missense"
    GNOMAD_API = "https://gnomad.broadinstitute.org/api"
    query = """
    query ($geneSymbol: String!) {
    gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
            variants(dataset: gnomad_r4) {
                variant_id
                pos
                ref
                alt
                consequence
                hgvsc
                hgvsp
            }
        }
    }
    """
    try:
        response = requests.post(GNOMAD_API, json={"query": query,"variables": {"geneSymbol": gene_symbol}})
        data = response.json()
        variants = data["data"]["gene"]["variants"]
        #gnomad_missense_list = [v for v in variants if "missense_variant" in v["consequence"]]
        return variants
    except requests.RequestException as e:
        print(f"Failed to fetch gnomAD missense variants for {gene_symbol}: {e}")
        return []
def fetch_all_clinvar_for_gene(gene, api_key=None):
    """
    Fetch all SNV records from ClinVar for a given gene symbol.
    Returns parsed XML root, or None on failure.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    search_params = {
        "db": "clinvar",
        "term": f"{gene}[gene] NOT \"copy number gain\"[Type] NOT \"copy number loss\"[Type]",
        "retmax": 5000,
        "retmode": "json",
        "usehistory": "y",
    }
    if api_key:
        search_params["api_key"] = api_key

    try:
        search_resp = requests.get(
            f"{base_url}/esearch.fcgi",
            params=search_params,
            timeout=30
        )
        search_resp.raise_for_status()
        esearch = search_resp.json().get("esearchresult", {})

        total = int(esearch.get("count", 0))
        if total == 0:
            return None, 0

        time.sleep(0.34 if not api_key else 0.11)

        fetch_params = {
            "db": "clinvar",
            "query_key": esearch.get("querykey"),
            "WebEnv": esearch.get("webenv"),
            "rettype": "vcv",
            "retmode": "xml",
            "retmax": 5000,
            "from_esearch": "yes",
        }
        if api_key:
            fetch_params["api_key"] = api_key

        fetch_resp = requests.get(
            f"{base_url}/efetch.fcgi",
            params=fetch_params,
            timeout=60
        )
        fetch_resp.raise_for_status()
        return ET.fromstring(fetch_resp.content), total

    except Exception as e:
        print(f"    Warning: ClinVar query failed for {gene}: {e}")
        return None, 0
def get_exon_intron_coor(transcript_id):
    url = f"https://rest.ensembl.org/lookup/id/{transcript_id}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    response = requests.get(url, headers=headers, params={"expand": 1, "utr": 1})
    response.raise_for_status()
    data = response.json()

    exons = sorted(data["Exon"], key=lambda e: e["start"])

    exon_map = [{"type": "exon", "start": e["start"], "end": e["end"]} for e in exons]

    # introns are the gaps between consecutive exons
    intron_map = [
        {"type": "intron", "start": exons[i]["end"] + 1, "end": exons[i+1]["start"] - 1}
        for i in range(len(exons) - 1)
    ]

    return sorted(exon_map + intron_map, key=lambda x: x["start"])
def get_cds_map(transcript_id):
    url = f"https://rest.ensembl.org/overlap/id/{transcript_id}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = {"feature": "CDS"}
    
    response = requests.get(url, headers=headers, params=params)
    response.raise_for_status()
    
    cds_regions = [{"start":i['start'],"end":i['end']} for i in response.json() if i ["Parent"]==transcript_id]
    cds_regions = sorted(cds_regions, key=lambda x: x["start"])

    return cds_regions
def adjust_display_map(display_map,intron_size):
    """
    Pass in a exon/intron map, with a desired intron_size to get back an adjusted exon/intron map with shrunken introns
    display_map: exon/intron map formated as list of [{'start':int, 'end':int},...]
    returns the display_map with shrunken introns
    """
    adj_display_map = []
    offset = display_map[0]['start']
    for region in display_map:
        if region['type'] == "exon":
            exon_width = region['end']-region['start']
            adj_display_map.append({'type':'exon','start':offset,'end':offset+exon_width})
            print(offset)
            offset += exon_width
        if region['type'] == "intron":
            adj_display_map.append({'type':'intron','start':offset,'end':offset+intron_size})
            offset += intron_size
    return adj_display_map
def adjust_trunc_variant_positions(trunc_var, display_map, intron_size):
    """
    Pass in the truncation variants list, display_map, and intron_size.
    This will return the truncation variants list with an adjusted_start and adjusted_end value
    """

    for region in display_map:
        region["width"] = region['end'] - region['start']

    adj_display_map = adjust_display_map(display_map,intron_size)

    adj_trunc_var_list = []
    for v in trunc_var:
        start, end = [int(i) for i in v['id'].split(":")[1].split("/")[0].split("-")]
        i = max([j for j in range(len(display_map)) if display_map[j]['start']<start])
        adj_start =adj_display_map[i]['start']+(start - display_map[i]['start'])
        adj_end =adj_display_map[i]['start']+(end - display_map[i]['start'])
        
        v['adj_start'] = adj_start
        v['adj_end']  = adj_end
    return trunc_var
def adjust_cds_map(cds_map, display_map, intron_size):

    adj_display_map = adjust_display_map(display_map,intron_size)
    offset = min([i['start'] for i in cds_map])
    for r in cds_map:
        start = r['start']
        end = r['end']
        i = max([j for j in range(len(display_map)) if display_map[j]['start']<=start])
        print(f"determined index: {i} brings cds_map[start]:{start} closest to display_map[start]{display_map[i]['start']}")
        adj_start = adj_display_map[i]['start']+(start - display_map[i]['start'])
        adj_end = adj_display_map[i]['start']+(start - display_map[i]['start'])+(end-start)

        r['adj_start'] = adj_start
        r['adj_end'] = adj_end
    return cds_map
def find_exon_number(pos, region_map):
    exons = [r for r in region_map if r["type"] == "exon"]
    for i, exon in enumerate(exons, start=1):
        if exon["start"] <= pos <= exon["end"]:
            return i
    return None  # in intron or not found
def plot_exon_map_plotly(display_map, clinvar_truncating_variants = None, gnomad_truncating_variants = None, cds_map = None, intron_size = None, box_height=1,gene_name = None, transcript_id = None):
    """
    display_map: output of build_display_map()
    truncating_variants: list of dicts with keys 'display_x' and any variant data to show on hover
    """
    fig = go.Figure()

    # Draw exons and introns
    x_exon, y_exon = [], []

    if intron_size is None:
        x_exon, y_exon = [], []
        x_intron, y_intron = [], []

        for region in display_map:
            if region["type"] == "exon":
                x_exon += [region["start"], region["start"], region["end"], region["end"], None]
                y_exon += [0, box_height, box_height, 0, None]
        # All introns as one line trace
            if region["type"] == "intron":
                mid = box_height / 2
                x_intron += [region["start"], region["end"], None]
                y_intron += [mid, mid, None]

    if not intron_size is None:
        adj_display_map = adjust_display_map(display_map,intron_size)
        x_exon, y_exon = [], []
        x_intron, y_intron = [], []

        for region in adj_display_map:
            if region['type'] == 'exon':
                x_exon += [region["start"], region["start"], region["end"], region["end"], None]
                y_exon += [0, box_height, box_height, 0, None]
            if region["type"] == "intron":
                mid = box_height / 2
                x_intron += [region["start"], region["end"], None]
                y_intron += [mid, mid, None]
    
    x_cds, y_cds = [], []
    if not cds_map is None:
        if intron_size is None:
            for c in cds_map:
                x_cds += [c['start'], c['start'],c['end'],c['end'],None]
                y_cds += [0, box_height, box_height, 0, None]
        if not intron_size is None:
            adj_cds_map = adjust_cds_map(cds_map,display_map,intron_size)
            for c in adj_cds_map:
                x_cds += [c['adj_start'], c['adj_start'],c['adj_end'],c['adj_end'],None]
                y_cds += [0, box_height, box_height, 0, None]
    fig.add_trace(go.Scatter(
        x=x_exon, y=y_exon,
        fill="toself", #fillcolor="steelblue",
        fillcolor = "purple",
        line=dict(color="purple"),
        mode="lines", showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=x_cds, y=y_cds,
        fill="toself", #fillcolor="steelblue",
        fillcolor = "steelblue",
        line=dict(color="steelblue"),
        mode="lines", showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=x_intron, y=y_intron,
        mode="lines", line=dict(color="black", width=1),
        showlegend=False, hoverinfo="skip"
    ))
    
    #================================================
    #adding the truncating variant marks
    #================================================
    #selects only the transcript_consequences for the transcript_id VEP annotations

    #adj_trunc_var = adjust_trunc_variant_positions(truncating_variants,display_map,intron_size)
    if not clinvar_truncating_variants is None:
        #selects only the transcript_consequences for the transcript_id VEP annotations
        clinvar_truncating_consequences = [item for sublist in clinvar_truncating_variants for item in sublist['transcript_consequences'] if item['transcript_id'] == transcript_id]

        for j in range(len(clinvar_truncating_variants)):
            clinvar_truncating_variants[j]['mane'] = clinvar_truncating_consequences[j]
        if intron_size is None:
            x_coor=[int(v["id"].split(":")[1].split("/")[0].split("-")[0]) for v in clinvar_truncating_variants]
            y_coor=[box_height/2 +random.uniform(0.3,0.8) for _ in clinvar_truncating_variants]

        if not intron_size is None:
            
            adj_trunc_var = adjust_trunc_variant_positions(clinvar_truncating_variants,display_map,intron_size)
            x_coor=[v['adj_start'] for v in adj_trunc_var]
            y_coor=[box_height/2 +random.uniform(0.3,0.8) for _ in adj_trunc_var]
        fig.add_trace(go.Scatter(
            x=x_coor, y=y_coor,
            yaxis ="y2", #puts variant traces on secondary y axis that isn't visible in slider
            mode="markers",
            marker=dict(color='red',size=8),
            #consequence_terms
            #amino_acids, protein_start, protein_end
            #impact
            customdata=[
                [v['id'], 
                f"p.{v['mane']['amino_acids']}-{v['mane']['protein_start']}-{v['mane']['protein_end']}", 
                v['mane']["consequence_terms"],
                v['mane']["impact"]]
                for v in clinvar_truncating_variants
            ],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "prot_loc: %{customdata[1]}<br>"
                "consequence: %{customdata[2]}<br>"
                "impact: %{customdata[3]}<br>"
                "<extra></extra>"
            ),
            name="ClinVar truncating variants"
        ))
    #print(adj_display_map)
    #================================================
    #adding the truncating variant marks from gnomad!!
    #================================================
    if not gnomad_truncating_variants is None:
        #selects only the transcript_consequences for the transcript_id VEP annotations
        gnomad_truncating_consequences = [item for sublist in gnomad_truncating_variants for item in sublist['transcript_consequences'] if item['transcript_id'] == transcript_id]
        for j in range(len(gnomad_truncating_variants)):
            gnomad_truncating_variants[j]['mane'] = gnomad_truncating_consequences[j]
        if intron_size is None:
            x_coor=[int(v["id"].split(":")[1].split("/")[0].split("-")[0]) for v in gnomad_truncating_variants]
            y_coor=[box_height/2 +random.uniform(-0.2,0.3) for _ in gnomad_truncating_variants]

        if not intron_size is None:
            
            adj_trunc_var = adjust_trunc_variant_positions(gnomad_truncating_variants,display_map,intron_size)
            x_coor=[v['adj_start'] for v in adj_trunc_var]
            y_coor=[box_height/2 +random.uniform(-0.2,0.3) for _ in adj_trunc_var]
        fig.add_trace(go.Scatter(
            x=x_coor, y=y_coor,
            yaxis ="y2", #puts variant traces on secondary y axis that isn't visible in slider
            mode="markers",
            marker=dict(color='blue',size=8),
            #consequence_terms
            #amino_acids, protein_start, protein_end
            #impact
            customdata=[
                [v['id'], 
                f"p.{v['mane']['amino_acids']}-{v['mane']['protein_start']}-{v['mane']['protein_end']}", 
                v['mane']["consequence_terms"],
                v['mane']["impact"]]
                for v in gnomad_truncating_variants
            ],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "prot_loc: %{customdata[1]}<br>"
                "consequence: %{customdata[2]}<br>"
                "impact: %{customdata[3]}<br>"
                "<extra></extra>"
            ),
            name="GnomAD truncating variants"
        ))



    x_min = (min(x_cds[::5]))-400
    x_max = (max(x_cds[2::5]))+400
    strand = list(set([v['strand'] for v in gnomad_truncating_variants]))[0]

    ########Layout configuration!!!
    fig.update_layout(
    
        xaxis=dict(title="Position (bp from start of chromosome)", showgrid=False, range=[x_min, x_max]),
        yaxis=dict(visible=False, range=[-0.5, box_height + 0.5]),
        yaxis2=dict(
            overlaying="y",
            visible=False,
            range=[-0.5, box_height + 0.5]  # match primary y range
        ),
        height=400,
        #width=1400,
        plot_bgcolor="white"
    )
    fig.update_layout(
        xaxis=dict(
            rangeslider=dict(visible=True),   # scrollable minimap below chart
            type="linear"
        )
    )
    fig.update_layout(
        xaxis=dict(showline=True, linewidth=2, linecolor="black", mirror=True),
        yaxis=dict(showline=True, linewidth=2, linecolor="black", mirror=True)
    )
    fig.update_layout(yaxis=dict(range=[-0.5, box_height + 3]))
         

    fig.update_layout(xaxis=dict(title="Position (bp from start of chromosome)", showgrid=False))
    fig.update_layout(
    xaxis=dict(
        autorange="reversed" if strand == -1 else True
    ),
    xaxis2=dict(  # minimap too
        autorange="reversed" if strand == -1 else True
    )
    )
    fig.update_layout(dragmode="pan")
    fig.update_layout(
        xaxis=dict(
            autorange=False,
            range=[x_max, x_min] if strand == -1 else [x_min, x_max]),
        xaxis2=dict(
            autorange=False,
            range=[x_max, x_min] if strand == -1 else [x_min, x_max]))  
     
    title=dict(text=f"{gene_name}: Early truncations: {transcript_id}",
            font=dict(size=24),
            x=0.5,          # center horizontally
            xanchor="center")

    return fig
def is_missense(variation):
    for mc in variation.findall(".//MolecularConsequence"):
        so_term  = mc.get("Type", "").lower()
        function = mc.get("Function", "").lower()
        if "missense" in so_term or "missense" in function:
            return True
    # HGVS protein fallback
    for hgvs in variation.findall(".//HGVSExpression[@Type='hgvs, protein']"):
        text = hgvs.text or ""
        if text and not any(x in text for x in ["Ter", "*", "fs", "del", "ins", "="]):
            return True
    return False
#==================================================================================
#for creating the gene dictionaries of missense variants
#==================================================================================
def create_individual_missense_gene_dict(gene):
    #get transcripts
    #canonical_transcript = get_uniprot_canonical_transcript(gene)
    MANE_transcript = get_mane_transcript(gene)
    
    
    #get sequences from transcripts
    prot_seqs = get_protein_sequences(gene, MANE_transcript)
    
    #get the binding/disorder scores
    binding_scores = {}
    disorder_scores = {}

    with open("postsynaptic_genes_aiupred_score_list.json", "r") as f:
        score_dict = json.load(f)
    score_dict = {k: v for d in score_dict for k, v in d.items()}
    binding_arr = score_dict[gene][transcript_id]['binding_scores']
    disorder_arr = score_dict[gene][transcript_id]['disorder_scores']
        
    #Get the ClinVar data
    root, tmp_variants = fetch_snvs_for_gene(gene, api_key=NCBI_API_KEY)
    if root is None:
        print(f"No ClinVar records found for {gene}. Skipping variant parsing.")
        #skipped_genes.append(gene)
        clinvar_missense_variants = []
        clinvar_density_curve = []  # No variants, so density is zero
    else:
        variations = root.findall(".//VariationArchive")
        clinvar_missense_variants = [
            my_parse_variants(v) for v in variations if is_missense(v)
        ]
    #Get the GnomAD data
    gnomad_missense_variants = fetch_gnomad_missense_variants(gene)  

    #Get the list of clinvar missense protein mutations by transcript
    clinvar_VEP_annotations = get_ensemblVEP_aa_consequences([i['variant_id'] for i in clinvar_missense_variants])
    clinvar_missense_prot_consq = assemble_ensemblVEP_aa_consequences(clinvar_VEP_annotations, MANE_transcript)

    #Get the list of gnomad missense protein mutations by transcript
    gnomad_VEP_annotations = get_ensemblVEP_aa_consequences([i['variant_id'] for i in gnomad_missense_variants])
    gnomad_missense_prot_consq = assemble_ensemblVEP_aa_consequences(gnomad_VEP_annotations, MANE_transcript)
        
    #Create clinvar density curve
    clinvar_density_curves = {}
    gnomad_density_curves = {}
    for transc_id, seq in prot_seqs.items():
        indexed_list = assemble_indexed_residue_list(clinvar_missense_prot_consq[transc_id],len(seq))
        clinvar_density_curves[transc_id] = convolve_v3(indexed_list, len(seq))
        indexed_list = assemble_indexed_residue_list(gnomad_missense_prot_consq[transc_id],len(seq))
        gnomad_density_curves[transc_id] = convolve_v3(indexed_list, len(seq))

    mane_clinvar_VEP = []
    for i in clinvar_VEP_annotations:
        for j in i['transcript_consequences']:
            if j['transcript_id']==MANE_transcript:
                mane_clinvar_VEP.append(j)
    mane_gnomad_VEP = []
    for i in gnomad_VEP_annotations:
        for j in i['transcript_consequences']:
            if j['transcript_id']==MANE_transcript:
                mane_gnomad_VEP.append(j)
    
    tmp_gene_dict = {}
    for transc_id in prot_seqs:
        tmp_gene_dict[transc_id] = {
            "sequence" :                    prot_seqs[transc_id],
            "length" :                      len(prot_seqs[transc_id]),
            "clinvar_missense_variants" :   clinvar_missense_prot_consq[transc_id],
            "gnomad_missense_variants" :    gnomad_missense_prot_consq[transc_id],
            "clinvar_density_curve" :       clinvar_density_curves[transc_id],
            "gnomad_density_curve" :        gnomad_density_curves[transc_id],
            "disorder_scores" :             disorder_arr,
            "binding_scores":               binding_arr
        }
    tmp_gene_dict["MANE_clinvar_VEP"] = mane_clinvar_VEP
    tmp_gene_dict["MANE_gnomad_VEP"] = mane_gnomad_VEP
    tmp_gene_dict["clinvar_VEP_annotations"] = clinvar_VEP_annotations
    tmp_gene_dict["gnomad_VEP_annotations"] = gnomad_VEP_annotations
    return tmp_gene_dict
def get_uniprot_canonical_transcript(gene_symbol):
    # Search for the reviewed (Swiss-Prot) human entry
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "fields": "xref_ensembl,xref_refseq,cc_alternative_products"
    }
    response = requests.get(url, params=params)
    response.raise_for_status()
    results = response.json()["results"]
    
    if not results:
        return None

    entry = results[0]
    transcripts = []

    canonical_isoform = [isoform['isoformIds'][0] for isoform in entry['comments'][0]['isoforms'] if isoform['isoformSequenceStatus'] == "Displayed"][0] #Get the protein isoform ID from Uniprot that has a Displayed status
    canonical_transcript = [t['id'] for t in entry["uniProtKBCrossReferences"] if 'isoformId' in t and t['isoformId'] == canonical_isoform]#using the ID of the canonical isoform, I check which transcripts transcribe that ID and return it as a list
    canonical_transcript = [i.split(".")[0] for i in canonical_transcript if i.startswith("ENST")]#find the transcripts that start with ENST and return them dropping anything after the "." (which indicates the version number)
    
    return canonical_transcript[0]






    """
    display_map: output of build_display_map()
    truncating_variants: list of dicts with keys 'display_x' and any variant data to show on hover
    """
    fig = go.Figure()

    # Draw exons and introns
    x_exon, y_exon = [], []

    if intron_size is None:
        x_exon, y_exon = [], []
        x_intron, y_intron = [], []

        for region in display_map:
            if region["type"] == "exon":
                x_exon += [region["start"], region["start"], region["end"], region["end"], None]
                y_exon += [0, box_height, box_height, 0, None]
        # All introns as one line trace
            if region["type"] == "intron":
                mid = box_height / 2
                x_intron += [region["start"], region["end"], None]
                y_intron += [mid, mid, None]

    if not intron_size is None:
        adj_display_map = adjust_display_map(display_map,intron_size)
        x_exon, y_exon = [], []
        x_intron, y_intron = [], []

        for region in adj_display_map:
            if region['type'] == 'exon':
                x_exon += [region["start"], region["start"], region["end"], region["end"], None]
                y_exon += [0, box_height, box_height, 0, None]
            if region["type"] == "intron":
                mid = box_height / 2
                x_intron += [region["start"], region["end"], None]
                y_intron += [mid, mid, None]
    
    x_cds, y_cds = [], []
    if not cds_map is None:
        if intron_size is None:
            for c in cds_map:
                x_cds += [c['start'], c['start'],c['end'],c['end'],None]
                y_cds += [0, box_height, box_height, 0, None]
        if not intron_size is None:
            adj_cds_map = adjust_cds_map(cds_map,display_map,intron_size)
            for c in adj_cds_map:
                x_cds += [c['adj_start'], c['adj_start'],c['adj_end'],c['adj_end'],None]
                y_cds += [0, box_height, box_height, 0, None]
    fig.add_trace(go.Scatter(
        x=x_exon, y=y_exon,
        fill="toself", #fillcolor="steelblue",
        fillcolor = "purple",
        line=dict(color="purple"),
        mode="lines", showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=x_cds, y=y_cds,
        fill="toself", #fillcolor="steelblue",
        fillcolor = "steelblue",
        line=dict(color="steelblue"),
        mode="lines", showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=x_intron, y=y_intron,
        mode="lines", line=dict(color="black", width=1),
        showlegend=False, hoverinfo="skip"
    ))
    
    #================================================
    #adding the truncating variant marks
    #================================================
    #selects only the transcript_consequences for the transcript_id VEP annotations

    #adj_trunc_var = adjust_trunc_variant_positions(truncating_variants,display_map,intron_size)
    if not clinvar_truncating_variants is None:
        #selects only the transcript_consequences for the transcript_id VEP annotations
        clinvar_truncating_consequences = [item for sublist in clinvar_truncating_variants for item in sublist['transcript_consequences'] if item['transcript_id'] == transcript_id]

        for j in range(len(clinvar_truncating_variants)):
            clinvar_truncating_variants[j]['mane'] = clinvar_truncating_consequences[j]
        if intron_size is None:
            x_coor=[int(v["id"].split(":")[1].split("/")[0].split("-")[0]) for v in clinvar_truncating_variants]
            y_coor=[box_height/2 +random.uniform(0.3,0.8) for _ in clinvar_truncating_variants]

        if not intron_size is None:
            
            adj_trunc_var = adjust_trunc_variant_positions(clinvar_truncating_variants,display_map,intron_size)
            x_coor=[v['adj_start'] for v in adj_trunc_var]
            y_coor=[box_height/2 +random.uniform(0.3,0.8) for _ in adj_trunc_var]
        fig.add_trace(go.Scatter(
            x=x_coor, y=y_coor,
            yaxis ="y2", #puts variant traces on secondary y axis that isn't visible in slider
            mode="markers",
            marker=dict(color='red',size=8),
            #consequence_terms
            #amino_acids, protein_start, protein_end
            #impact
            customdata=[
                [v['id'], 
                f"p.{v['mane']['amino_acids']}-{v['mane']['protein_start']}-{v['mane']['protein_end']}", 
                v['mane']["consequence_terms"],
                v['mane']["impact"]]
                for v in clinvar_truncating_variants
            ],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "prot_loc: %{customdata[1]}<br>"
                "consequence: %{customdata[2]}<br>"
                "impact: %{customdata[3]}<br>"
                "<extra></extra>"
            ),
            name="ClinVar truncating variants"
        ))
    #print(adj_display_map)
    #================================================
    #adding the truncating variant marks from gnomad!!
    #================================================
    if not gnomad_truncating_variants is None:
        #selects only the transcript_consequences for the transcript_id VEP annotations
        gnomad_truncating_consequences = [item for sublist in gnomad_truncating_variants for item in sublist['transcript_consequences'] if item['transcript_id'] == transcript_id]
        for j in range(len(gnomad_truncating_variants)):
            gnomad_truncating_variants[j]['mane'] = gnomad_truncating_consequences[j]
        if intron_size is None:
            x_coor=[int(v["id"].split(":")[1].split("/")[0].split("-")[0]) for v in gnomad_truncating_variants]
            y_coor=[box_height/2 +random.uniform(-0.2,0.3) for _ in gnomad_truncating_variants]

        if not intron_size is None:
            
            adj_trunc_var = adjust_trunc_variant_positions(gnomad_truncating_variants,display_map,intron_size)
            x_coor=[v['adj_start'] for v in adj_trunc_var]
            y_coor=[box_height/2 +random.uniform(-0.2,0.3) for _ in adj_trunc_var]
        fig.add_trace(go.Scatter(
            x=x_coor, y=y_coor,
            yaxis ="y2", #puts variant traces on secondary y axis that isn't visible in slider
            mode="markers",
            marker=dict(color='blue',size=8),
            #consequence_terms
            #amino_acids, protein_start, protein_end
            #impact
            customdata=[
                [v['id'], 
                f"p.{v['mane']['amino_acids']}-{v['mane']['protein_start']}-{v['mane']['protein_end']}", 
                v['mane']["consequence_terms"],
                v['mane']["impact"]]
                for v in gnomad_truncating_variants
            ],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "prot_loc: %{customdata[1]}<br>"
                "consequence: %{customdata[2]}<br>"
                "impact: %{customdata[3]}<br>"
                "<extra></extra>"
            ),
            name="GnomAD truncating variants"
        ))



    x_min = (min(x_cds[::5]))-400
    x_max = (max(x_cds[2::5]))+400
    strand = list(set([v['strand'] for v in gnomad_truncating_variants]))[0]

    ########Layout configuration!!!
    fig.update_layout(
    
        xaxis=dict(title="Position (bp from start of chromosome)", showgrid=False, range=[x_min, x_max]),
        yaxis=dict(visible=False, range=[-0.5, box_height + 0.5]),
        yaxis2=dict(
            overlaying="y",
            visible=False,
            range=[-0.5, box_height + 0.5]  # match primary y range
        ),
        height=400,
        #width=1400,
        plot_bgcolor="white"
    )
    fig.update_layout(
        xaxis=dict(
            rangeslider=dict(visible=True),   # scrollable minimap below chart
            type="linear"
        )
    )
    fig.update_layout(
        xaxis=dict(showline=True, linewidth=2, linecolor="black", mirror=True),
        yaxis=dict(showline=True, linewidth=2, linecolor="black", mirror=True)
    )
    fig.update_layout(yaxis=dict(range=[-0.5, box_height + 3]))
         

    fig.update_layout(xaxis=dict(title="Position (bp from start of chromosome)", showgrid=False))
    fig.update_layout(
    xaxis=dict(
        autorange="reversed" if strand == -1 else True
    ),
    xaxis2=dict(  # minimap too
        autorange="reversed" if strand == -1 else True
    )
    )


    fig.update_layout(
        xaxis=dict(
            autorange=False,
            range=[x_max, x_min] if strand == -1 else [x_min, x_max]),
        xaxis2=dict(
            autorange=False,
            range=[x_max, x_min] if strand == -1 else [x_min, x_max]))  
     
    title=dict(text=f"{gene_name}: Early truncations: {transcript_id}",
            font=dict(size=24),
            x=0.5,          # center horizontally
            xanchor="center")

    return fig
def get_protein_sequences(gene_name, *transcript_ids):
    """
    Get amino acid sequences for a gene across multiple transcripts.

    Args:
        gene_name (str): Gene symbol for reference.
        *transcript_ids (str): Any number of Ensembl transcript IDs (ENST...).

    Returns:
        dict: Maps each transcript_id to its amino acid sequence string.
    """
    url = "https://rest.ensembl.org/sequence/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    results = {}
    for transcript_id in transcript_ids:
        response = requests.get(
            f"{url}/{transcript_id}",
            headers=headers,
            params={"type": "protein"}
        )
        response.raise_for_status()
        results[transcript_id] = response.json()["seq"]
    
    return results
def get_aiupred_scores(sequence):
    #This function will take a protein sequence as input and return the aiupred binding score and disorder score for each residue in the sequence
    predictor = AIUPredMPS()
    binding_scores = predictor.predict_binding(sequence).tolist()
    disorder_scores = predictor.predict_disorder(sequence).tolist()
    return binding_scores, disorder_scores
def fetch_snvs_for_gene(gene, api_key=None):
    """
    Fetch all SNV records from ClinVar for a given gene symbol.
    Returns parsed XML root, or None on failure.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    search_params = {
        "db": "clinvar",
        "term": f"{gene}[gene] AND \"single nucleotide variant\"[Type]",
        "retmax": 5000,
        "retmode": "json",
        "usehistory": "y",
    }
    if api_key:
        search_params["api_key"] = api_key

    try:
        search_resp = requests.get(
            f"{base_url}/esearch.fcgi",
            params=search_params,
            timeout=30
        )
        search_resp.raise_for_status()
        esearch = search_resp.json().get("esearchresult", {})

        total = int(esearch.get("count", 0))
        if total == 0:
            return None, 0

        time.sleep(0.34 if not api_key else 0.11)

        fetch_params = {
            "db": "clinvar",
            "query_key": esearch.get("querykey"),
            "WebEnv": esearch.get("webenv"),
            "rettype": "vcv",
            "retmode": "xml",
            "retmax": 5000,
            "from_esearch": "yes",
        }
        if api_key:
            fetch_params["api_key"] = api_key

        fetch_resp = requests.get(
            f"{base_url}/efetch.fcgi",
            params=fetch_params,
            timeout=60
        )
        fetch_resp.raise_for_status()
        return ET.fromstring(fetch_resp.content), total

    except Exception as e:
        print(f"    Warning: ClinVar query failed for {gene}: {e}")
        return None, 0
def is_truncation(variation):
    truncating = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "nonsense"}

    # How to determine if an early truncation
    # HGVS protein fallback
    for hgvs in variation.findall(".//HGVSExpression[@Type='hgvs, protein']"):
        text = hgvs.text or ""
        if text and any(x in text for x in ["Ter", "*", "fs", "del", "ins"]):
            return True
    return False
def fetch_gnomad_missense_variants(gene_symbol):
    #This function will query gnomAD for missense variants in the given gene symbol and return a list of missense variants with their positions and allele frequencies
    #I can use the gnomAD API to fetch this information. The endpoint for fetching variants by gene is "https://gnomad.broadinstitute.org/api/variants/search?query=gene:{gene_symbol}&variantType=missense"
    GNOMAD_API = "https://gnomad.broadinstitute.org/api"
    query = """
    query ($geneSymbol: String!) {
    gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
            variants(dataset: gnomad_r4) {
                variant_id
                pos
                ref
                alt
                consequence
                hgvsc
                hgvsp
            }
        }
    }
    """
    try:
        response = requests.post(GNOMAD_API, json={"query": query,"variables": {"geneSymbol": gene_symbol}})
        data = response.json()
        variants = data["data"]["gene"]["variants"]
        gnomad_missense_list = [v for v in variants if "missense_variant" in v["consequence"]]
        return gnomad_missense_list
    except requests.RequestException as e:
        print(f"Failed to fetch gnomAD missense variants for {gene_symbol}: {e}")
        return []
def assemble_ensemblVEP_aa_consequences(total_ensemblVEP_consequences,*transcript_ids):
    results = []
    #I want a list of key:value pairs. The key is the transcript_ID and the pair is the missense variants
    for v in total_ensemblVEP_consequences: 
        tc = v['transcript_consequences']
        tc = [i for i in tc if i['transcript_id'] in transcript_ids]
        for t in tc:
            aa = t.get("amino_acids","N/A")
            if aa != "N/A":
                #print("processing amino acid: "+ str(aa))
                results.append({t['transcript_id']:str(aa[0])+str(t['protein_start'])+str(aa[-1])})
    
    tmp_results = {}
    for t in transcript_ids:
        tmp_results[t]  = ([r[t] for r in results if list(r.keys())[0] == t])
    results = tmp_results
    
    print("successfully fetched all protein consequences for the transcripts: "+str([i for i in transcript_ids]))
    return results
def assemble_indexed_residue_list(residues,prot_len):
    indexed_residues = [0]*prot_len
    residues_num = [int(i[1:-1]) for i in residues]
    print(residues_num)
    for i in range(prot_len):
        if i in residues_num:
            indexed_residues[i] = 1
    return indexed_residues  
def convolve_v3(rsd_index_list, protein_len, window_size=5):
    n = window_size*2
    p = 0.5
    r_values = list(range(n+1))
    dist = [binom.pmf(r,n,p) for r in r_values ]
    # plt.bar(r_values, dist)
    # plt.show()
    norm = [float(i)/max(dist) for i in dist]
    
    convolved_list = []
    for i in range(protein_len):
        tmp = 0
        
        offset_numerator = max([(0 if (i-protein_len+int(len(dist)/2)) <0 else (i-protein_len+int(len(dist)/2))), 
                             (0 if -1*(i-int(len(dist)/2)) < 0 else -1*(i-int(len(dist)/2)))])
        #print(offset_numerator)
        
        for j in range(len(dist)):
            #if i - int(len(dist)/2) >0:
            #checks if the j values are within the range of the i values I guess *** go back and check this
            if i-int(len(dist)/2)+j >= 0 and i-int(len(dist)/2)+j < protein_len:
                tmp += rsd_index_list[i-int(len(dist)/2)+j]*dist[j]
            
        convolved_list.append(tmp*(1+(offset_numerator/(len(dist)/2))))
    return convolved_list

def plot_missense_cluster_chart(gene_dict,gene):

    DISORDER_THRESHOLD = 0.5
    transcript_id = get_mane_transcript(gene)
    data = gene_dict[transcript_id]
    length = data["length"]

    # with open("postsynaptic_genes_aiupred_score_list.json", "r") as f:
    #     score_dict = json.load(f)
    # binding_arr = np.array(score_dict[0][gene][transcript_id]['binding_scores'])
    # disorder_arr = np.array(score_dict[0][gene][transcript_id]['disorder_scores'])
    binding_arr = data['binding_scores']
    disorder_arr = data['disorder_scores']

    box_height=1
    clinvar_missense = data['clinvar_missense_variants']
    gnomad_missense = data['gnomad_missense_variants']
    clinvar_density = data['clinvar_density_curve']
    gnomad_density = data['gnomad_density_curve']

    x = list(range(1, length + 1))
    y_diff_arr = np.array(np.subtract(gnomad_density, clinvar_density))

    #fig = make_subplots(rows=1, cols=1, row_heights=[1],
    #                    shared_xaxes=True, vertical_spacing=0.02)
    fig = go.Figure()

    y_diff_trace = go.Scatter(x=x, y=list(y_diff_arr),mode="lines", name="tolerance",line=dict(color="purple"),legendgroup='tolerance')
    fill_blue_trace = go.Scatter(x=x, y=[v if v > 0 else 0 for v in y_diff_arr],fill="tozeroy", fillcolor="rgba(0,0,255,0.4)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='tolerance')
    fill_red_trace =go.Scatter(x=x, y=[v if v < 0 else 0 for v in y_diff_arr],fill="tozeroy", fillcolor="rgba(255,0,0,0.4)",
        line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip", legendgroup='tolerance')

    fill_clinvar_density_trace =go.Scatter(x=x, y=clinvar_density,fill="tozeroy", name = 'ClinVar densities', fillcolor="rgba(255,0,0,0.2)",
        line=dict(color="red"),showlegend=True, hoverinfo="skip", legendgroup='ClinVar densities', visible = 'legendonly')

    fill_gnomad_density_trace =go.Scatter(x=x, y=gnomad_density,fill="tozeroy", name = 'GnomAD densities', fillcolor="rgba(0,0,255,0.2)",
        line=dict(color="blue"),showlegend=True, hoverinfo="skip", legendgroup='GnomAD densities', visible = 'legendonly')

    binding_score_trace = go.Scatter(x=x, y=list(binding_arr),mode="lines", name="AIUPred Binding",
        line=dict(color="green"),legendgroup='aiupred_binding', visible = 'legendonly')
    fill_binding_pink_trace = go.Scatter(x=x, y=[v if v >= DISORDER_THRESHOLD else 0 for v in binding_arr],fill="tozeroy", 
        fillcolor="rgba(255, 192, 203, 0.3)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='aiupred_binding',visible = 'legendonly')
    
    disorder_score_trace = go.Scatter(x=x, y=list(disorder_arr),mode="lines", name="AIUPred Disorder",
        line=dict(color="yellow"),legendgroup='aiupred_disorder', visible = 'legendonly')
    fill_disorder_green_trace = go.Scatter(x=x, y=[v if v >= DISORDER_THRESHOLD else 0 for v in disorder_arr],fill="tozeroy", 
        fillcolor="rgba(0,128,0,0.3)",line=dict(color="rgba(0,0,0,0)"),showlegend=False, hoverinfo="skip",legendgroup='aiupred_disorder',visible = 'legendonly')
 
    clinvar_missense_y = [-1*((i % 8) * (1 / 20) + 0.3) for i in range(len(clinvar_missense))]
    clinvar_missense_scatter = go.Scatter(x = [i[1:-1]for i in clinvar_missense],y = [-1*((i % 8) * (1 / 20) + 0.3) for i in range(len(clinvar_missense))],
            yaxis ="y2", #puts variant traces on secondary y axis that isn't visible in slider
            mode="markers",
            marker=dict(color="red",size=5),
            customdata=[v for v in clinvar_missense],
            hovertemplate=(
                "<b>%{customdata}</b><br>"
                "<extra></extra>"
            ),
            name="ClinVar")
    gnomad_missense_y = [-1*((i % 8) * (1 / 20) + 0.3) for i in range(len(gnomad_missense))]
    gnomad_missense_scatter = go.Scatter(
            x = [i[1:-1]for i in clinvar_missense],
            #y = [-1*(box_height/2 +random.uniform(-0.2,0.3)) for _ in clinvar_missense],
            y = [(i % 8) * (1 / 20) + 0.3 for i in range(len(gnomad_missense))],

            yaxis ="y2", #puts variant traces on secondary y axis that isn't visible in slider
            mode="markers",
            marker=dict(color="blue",size=5),
            customdata=[v for v in clinvar_missense],
            hovertemplate=(
                "<b>%{customdata}</b><br>"
                "<extra></extra>"
            ),
            name="GnomAD")


    fig.add_trace(y_diff_trace)
    fig.add_trace(fill_blue_trace)
    fig.add_trace(fill_red_trace)
    fig.add_trace(clinvar_missense_scatter)
    fig.add_trace(gnomad_missense_scatter)

    fig.add_trace(fill_clinvar_density_trace)
    fig.add_trace(fill_gnomad_density_trace)

    fig.add_trace(binding_score_trace)
    fig.add_trace(fill_binding_pink_trace)
    fig.add_trace(disorder_score_trace)
    fig.add_trace(fill_disorder_green_trace)

    fig.update_layout(
        yaxis=dict(range=[min(y_diff_arr)-0.3, max(y_diff_arr)+0.3], fixedrange=True),
        yaxis2=dict(overlaying="y", visible=False, fixedrange=True, range=[min(y_diff_arr)-0.1, max(y_diff_arr)+0.1])
    )
    fig.update_layout(dragmode="pan")
    # fig.update_layout(
    #     xaxis2=dict(
    #         rangeslider=dict(visible=True, thickness=0.2,range=[0,len(binding_arr)])
    # ))
    fig.update_xaxes(rangeslider_visible=True)
    fig.update_xaxes(
        range=[0, len(binding_arr)],
        autorange=False,
        rangeslider=dict(
            visible=True,
            range=[0, len(binding_arr)]  # lock the slider range too
        )
    )
    return fig

def get_last_coding_exon_number(transcript_id,region_map = None, cds_map = None):
    if region_map is None:
        region_map = get_exon_intron_coor(transcript_id)
    if cds_map is None:
        cds_regions = get_cds_map(transcript_id)

    exons = [r for r in region_map if r["type"] == "exon"]

    last_coding_exon = None
    for i, exon in enumerate(exons, start=1):
        for cds in cds_map:
            if cds["start"] <= exon["end"] and cds["end"] >= exon["start"]:
                last_coding_exon = i
                break
    return last_coding_exon

def highlight_by_database(val):
    if val == 'ClinVar':
            return "background-color: red"
    if val == 'GnomAD':
        return "background-color: blue"

def create_truncation_df(trunc_clinvar, trunc_gnomad, transcript_id):
    clinvar_df = []
    region_map = get_exon_intron_coor(transcript_id)
    cds_map = get_cds_map(transcript_id)
    for i in trunc_clinvar:
        for j in i['transcript_consequences']:
            if j['transcript_id'] == transcript_id:
                exon_num = find_exon_number(int(i['id'].split(':')[1].split('/')[0].split('-')[0]), region_map)
                tmp = {'id':i['id'],
                    'gene_start' : int(i['id'].split(':')[1].split('/')[0].split('-')[0]),
                    'gene_end': int(i['id'].split(':')[1].split('/')[0].split('-')[1]),
                    'codons':j['codons'], 
                    'exon':exon_num,
                    'protein start':j['protein_start'],
                    'amino_acids':j['amino_acids'],
                    'database':'ClinVar',
                    'cds_start':j['cds_start'],
                    'consequence_terms':j['consequence_terms'],
                    
                    }
                clinvar_df.append(tmp)
                #print(f"appending:{tmp}")
    clinvar_df = pd.DataFrame(clinvar_df)
    gnomad_df = []
    for i in trunc_gnomad:
        for j in i['transcript_consequences']:
            if j['transcript_id'] == transcript_id:
                exon_num = find_exon_number(int(i['id'].split(':')[1].split('/')[0].split('-')[0]), region_map)
                tmp = {'id':i['id'],
                    'gene_start' : int(i['id'].split(':')[1].split('/')[0].split('-')[0]),
                    'gene_end': int(i['id'].split(':')[1].split('/')[0].split('-')[1]),
                    'codons':j['codons'], 
                    'exon':exon_num,
                    'protein start':j['protein_start'],
                    'amino_acids':j['amino_acids'],
                    'database':'GnomAD',
                    'cds_start':j['cds_start'],
                    'consequence_terms':j['consequence_terms'],
                    
                    }
                gnomad_df.append(tmp)
                #print(f"appending:{tmp}")

    gnomad_df = pd.DataFrame(gnomad_df)
    df = pd.concat([clinvar_df, gnomad_df], ignore_index = True)
    return df

def create_missense_df(missense_dict,transcript_id):
    
    for i in missense_dict['MANE_clinvar_VEP']:
        i['database'] = "ClinVar"
    for i in missense_dict["MANE_gnomad_VEP"]:
        i['database'] = "GnomAD"
    df = []
    for i in missense_dict['MANE_clinvar_VEP']:
        if 'missense_variant' in i['consequence_terms']:
            tmp = {'position':i['protein_start'],
                   'variation': i['amino_acids'],
                   'cDNA Position': i['cdna_start'],
                   'codons': i['codons'],
                   'database' : i['database'],
                   'consequence': i['consequence_terms'],
                   'polyphen prediction': i['polyphen_prediction'],
                   'polyphen score': i['polyphen_score'],
                   'sift prediction': i['sift_prediction'],
                   'sift score': i['sift_score'],
                   'impact':i['impact']}
            df.append(tmp)

    for i in missense_dict['MANE_gnomad_VEP']:
        if 'missense_variant' in i['consequence_terms']:
            tmp = {'position':i['protein_start'],
                   'variation': i['amino_acids'],
                   'cDNA Position': i['cdna_start'],
                   'codons': i['codons'],
                   'database' : i['database'],
                   'consequence': i['consequence_terms'],
                   'polyphen prediction': i['polyphen_prediction'],
                   'polyphen score': i['polyphen_score'],
                   'sift prediction': i['sift_prediction'],
                   'sift score': i['sift_score'],
                   'impact':i['impact']}
            df.append(tmp)

    # clinvar_df = pd.DataFrame(missense_dict['MANE_clinvar_VEP'])
    # gnomad_df = pd.DataFrame(missense_dict['MANE_gnomad_VEP'])
    # df = pd.concat([clinvar_df,gnomad_df],ignore_index=True)
    df = pd.DataFrame(df)
    return df


st.set_page_config(layout="wide")  # must be the first st. call in your script
st.title("Roche Lab Variant Viewer 1.0.0-alpha.1(May, 2026")

# gene = st.text_input("Gene name")
# transcript_id = st.text_input("Transcript ID", "MANE")
# intron_size = st.text_input("Adjust Intron Size", "200")
# Input row at the top
col1, col2, col3 = st.columns(3)
with col1:
    gene = st.text_input("Gene name")
with col2:
    transcript_id = st.text_input("Transcript ID", "MANE")
with col3:
    intron_size = st.text_input("Adjust Intron Size", "200")


if st.button("Generate"):
    with st.status("Generating chart (this will take a minute or two while retrieving all variant information)...", expanded=True) as status:
            #Fetch and write the truncation variants
        if transcript_id == "MANE":
            st.write("Fetching MANE transcipt_id")
            transcript_id = get_mane_transcript(gene)

        # if os.path.exists(f"{gene}-{transcript_id}_clinvar_early_truncations.json"):
        #     st.write("Locally accessing clinvar_early_truncation variants")
        #     with open(f"{gene}-{transcript_id}_clinvar_early_truncations.json", "r") as f:
        #         my_dict = json.load(f)
        #         trunc_clinvar = my_dict

        # if os.path.exists(f"{gene}-{transcript_id}_gnomad_early_truncations.json"):
        #     st.write("Locally accessing clinvar_early_truncation variants")
        #     with open(f"{gene}-{transcript_id}_gnomad_early_truncations.json", "r") as f:
        #         my_dict = json.load(f)
        #         trunc_gnomad = my_dict

        # if not (os.path.exists(f"{gene}-{transcript_id}_gnomad_early_truncations.json") or os.path.exists(f"{gene}-{transcript_id}_clinvar_early_truncations.json")):
        #     st.write("Fetching variants from clinvar and gnomad, filtering for truncations, and collating...")
        #     trunc_clinvar, trunc_gnomad = collate_truncation_variants(gene)

        st.write("Fetching variants from clinvar and gnomad, filtering for truncations, and collating...")
        trunc_clinvar, trunc_gnomad = collate_truncation_variants(gene)
        st.write("Successfully fetched and collated all truncation variants...")
            # with open(f"{gene}-{transcript_id}_clinvar_early_truncations.json", "w") as f:
            #     json.dump(trunc_clinvar, f, indent=2)
            # with open(f"{gene}-{transcript_id}_gnomad_early_truncations.json", "w") as f:
            #     json.dump(trunc_gnomad, f, indent=2)

        # if os.path.exists(f"{gene}-missense.json"):
        #     st.write("Locally accessing clinvar/gnomad missense variants")
        #     with open(f"{gene}-missense.json", "r") as f:
        #         my_dict = json.load(f)
        #         missense_dict = my_dict

        # if not (os.path.exists(f"{gene}-missense.json")):
        #     st.write("Fetching missense variants from clinvar and gnomad, filtering for truncations, and collating...")
        #     missense_dict = create_individual_missense_gene_dict(gene)
        #     with open(f"{gene}-missense.json", "w") as f:
        #         json.dump(missense_dict, f, indent=2)
        st.write("Fetching missense variants from clinvar and gnomad, filtering, and collating...")
        missense_dict = create_individual_missense_gene_dict(gene)


        st.write("Fetching exon/intron map...")
        region_map = get_exon_intron_coor(transcript_id)
        st.write("Fetching coding region coordinates...")
        cds_regions = get_cds_map(transcript_id)


        try:
            intron_size = int(intron_size)
            st.write(f"Adjusting intron_size to {intron_size}...")
        except ValueError:
            intron_size = None

        st.write(f"Creating chart...")
        fig = plot_exon_map_plotly(get_exon_intron_coor(transcript_id), trunc_clinvar, trunc_gnomad, intron_size = intron_size, cds_map = get_cds_map(transcript_id),gene_name = gene, transcript_id = transcript_id)
        fig_missense = plot_missense_cluster_chart(missense_dict,gene)
        fig.update_layout(height = 600)
        status.update(label="Done!", state="complete", expanded=False)

    last_exon = get_last_coding_exon_number(transcript_id, region_map, cds_regions)
    def highlight_last_exon_cell(val):
        return "background-color: yellow" if val == last_exon else""
    tab1, tab2 = st.tabs(["Truncations", "Missense"])
    tab1.plotly_chart(fig, use_container_width=True,
                    config={"scrollZoom": True})
    
    tab2.plotly_chart(fig_missense, use_container_width=True,
                    config={"scrollZoom": True})
    
    df_trunc = create_truncation_df(trunc_clinvar, trunc_gnomad, transcript_id)
    tab1.dataframe(df_trunc.style.map(highlight_last_exon_cell, subset=['exon']).map(highlight_by_database, subset=['database']))
    
    df_missense = create_missense_df(missense_dict,transcript_id)
    tab2.dataframe(df_missense.style.map(highlight_by_database, subset=['database']))

    
    #st.dataframe(df, use_container_width=True, height = 300)

    #==================================================
    #Add the protein structure element to my page
    #Next step: take the PDB structure and create the file pseudo-colored w the tolerance values
    #==================================================
    # import streamlit.components.v1 as components

    # mmdb_id = "8GS3"  # swap this for your protein's MMDB ID
    
    # icn3d_html = f"""
    # <iframe 
    #     allow="xr-spatial-tracking *" 
    #     src="https://www.ncbi.nlm.nih.gov/Structure/icn3d/?mmdbid={mmdb_id}&width=600&height=500&showcommand=0&mobilemenu=1&showtitle=0" 
    #     width="620" 
    #     height="520" 
    #     style="border:none">
    # </iframe>
    # """

    # components.html(icn3d_html, height=540)
