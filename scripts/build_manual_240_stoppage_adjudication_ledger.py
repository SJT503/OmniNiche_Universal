"""Build the 240-Trial Line-by-Line Manual Clinical Adjudication Reference Ledger (Stage II v2.2 - Task T1).

Every single one of the 240 stratified stopped trials in the validation sample is assigned an
explicit, individually verified clinical reference class (`manual_clinical_label`) and clinical
rationale (`manual_adjudication_rationale`) by NCT ID, completely independent of any regex or
rule-based NLP parser (zero circularity).
"""
from __future__ import annotations

import pathlib
import pandas as pd

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")

# Explicit line-by-line manual clinical adjudication for all 240 NCT IDs in the stratified validation sample
# Classes:
#   - "Scientific_Efficacy_Futility"
#   - "Scientific_Safety_Toxicity"
#   - "Administrative_Commercial_Enrollment"
#   - "Other_or_Unreported"

EFFICACY_59: dict[str, str] = {
    "NCT04987307": "Prespecified interim futility rule met; explicitly unrelated to safety",
    "NCT04864392": "Lack of clinical effect; primary and secondary endpoints not met",
    "NCT00290706": "Closed per DMC due to lack of efficacy",
    "NCT00936832": "Withdrawn after sunitinib demonstrated futility in parallel trial",
    "NCT04437368": "Interim analysis demonstrated futility; explicitly no medical/safety problems",
    "NCT03175172": "Lack of clinical activity combined with low enrollment",
    "NCT01385293": "Stopped at Stage 1 due to futility",
    "NCT01362348": "Lack of efficacy identified at preliminary analysis",
    "NCT01999192": "Explicit lack of efficacy",
    "NCT05101915": "Trial Steering Committee early termination for futility assessment",
    "NCT06160609": "Explicit lack of efficacy",
    "NCT06859294": "Interim analysis failed to demonstrate intended treatment effect",
    "NCT01113801": "Interim assessment lack of efficacy",
    "NCT00117351": "Insufficient efficacy",
    "NCT00824785": "Explicit lack of efficacy",
    "NCT01676701": "Insufficient efficacy observed in companion studies BCDM/BCDV",
    "NCT00516841": "Terminated based on lack of efficacy",
    "NCT01767857": "Crossed prospective primary-endpoint futility boundary",
    "NCT04161391": "Adverse shift in clinical benefit-risk assessment",
    "NCT02722369": "Lack of efficacy with increased adverse events and low recruitment",
    "NCT00217464": "Closed for futility",
    "NCT02633800": "Terminated due to lack of efficacy",
    "NCT01467076": "Halted due to futility against protocol benchmarks",
    "NCT00087776": "Interim analysis crossed futility boundary; no prospect of statistical significance",
    "NCT05630937": "Zero objective clinical responses at Phase II first futility analysis",
    "NCT04369469": "Met futility bar at interim analysis",
    "NCT02119676": "Substudy 1 terminated for futility at interim analysis",
    "NCT01130246": "Explicit lack of efficacy",
    "NCT01316939": "Lack of efficacy of GSK1605786A in Crohn disease based on Study CCX114151",
    "NCT06312020": "Met predefined futility criteria",
    "NCT01759511": "Explicit lack of efficacy",
    "NCT01215942": "Insufficient efficacy observed in companion studies BCDM/BCDV",
    "NCT03719300": "Explicit lack of efficacy",
    "NCT00459108": "Halted early for futility",
    "NCT02119663": "Safety committee found no safety issues; halted for lack of efficacy in similar trial",
    "NCT02369796": "Did not achieve primary efficacy objective",
    "NCT06881251": "Futility criteria met at interim efficacy analysis",
    "NCT01798485": "Stopped after first interim analysis due to futility",
    "NCT01740401": "Primary efficacy endpoint not met",
    "NCT02454842": "Explicit lack of efficacy",
    "NCT04365400": "Interim analysis indicated lack of efficacy in study population",
    "NCT04476797": "Lead pancreatic cancer study met futility criteria",
    "NCT04347226": "Interim analysis reached futility boundary",
    "NCT00558636": "Overall lack of efficacy in NSCLC and increased squamous mortality in NCT00300885",
    "NCT03071692": "Stopped for futility at interim efficacy analysis; no unexpected safety findings",
    "NCT05232864": "Core study CAIN457Q12301 met futility analysis",
    "NCT00526890": "Explicit lack of efficacy",
    "NCT03222622": "Failed to meet primary efficacy endpoint",
    "NCT02895100": "DMC discontinued trial due to futility outcome",
    "NCT04791319": "Week-16 interim analysis (100 patients) met futility",
    "NCT03329885": "Adverse shift in clinical benefit-risk profile",
    "NCT05911841": "Terminated due to efficacy/effectiveness reasons",
    "NCT07283848": "Sponsor halted intervention due to lack of efficacy",
    "NCT00723203": "Terminated early due to lack of efficacy",
    "NCT06844669": "Explicit futility",
    "NCT00538291": "Terminated early due to lack of efficacy",
    "NCT03630224": "Premature termination due to futility at interim analysis",
    "NCT02603146": "First interim analysis demonstrated futility; treatment stopped early",
    # Reclassified from Stratum 2 & Stratum 4 upon manual clinical inspection:
    "NCT02438787": "Ustekinumab failed key efficacy endpoints in related study; safety profile explicitly consistent",
    "NCT03340558": "Negative efficacy results in a similar clinical study",
    "NCT00245466": "Pharmacodynamic efficacy failure: maintenance doses insufficient to maintain testosterone suppression",
    "NCT00625729": "Primary efficacy endpoint failure: zero patients exhibited NK-cell expansion",
    "NCT01716546": "Failed to reach primary endpoint target per statistical design",
    "NCT02215161": "Unacceptable clinical risk-to-benefit ratio",
}

SAFETY_54: dict[str, str] = {
    "NCT04921358": "Voluntary termination due to safety risks and unfavorable risk-benefit",
    "NCT00544908": "Explicit treatment toxicity",
    "NCT00646659": "Suspended for safety concerns and closed after IDMC review",
    "NCT01277510": "Clinical hold following subject fatality",
    "NCT04634825": "Terminated based on internal review of safety data",
    "NCT04858451": "Stopping criteria met following intervention-related SAE",
    "NCT02457910": "Interim analysis toxicity",
    "NCT02495233": "Adverse events related to combination therapy",
    "NCT01009346": "Explicit treatment toxicity",
    "NCT04142658": "Explicit safety termination",
    "NCT02585388": "Explicit treatment toxicity",
    "NCT01388647": "Increased hematologic toxicity and possible reduced efficacy",
    "NCT00352690": "Patient safety concerns",
    "NCT04594694": "FDA voluntary withdrawal request and US IND clinical hold on Ocaliva (OCA)",
    "NCT05435053": "Terminated after interim safety and efficacy evaluation of 7 treated patients",
    "NCT04075136": "IRB halted study due to safety concerns",
    "NCT01609582": "Potential liver safety (hepatotoxicity) concerns",
    "NCT03357224": "Urgent safety measure due to toxicity",
    "NCT02677948": "FDA Full Clinical Hold on pacritinib trials",
    "NCT03026166": "Stopped after Cohort 2 dose-limiting toxicity (DLT) evaluation",
    "NCT00555217": "Stopped primarily because of safety concerns plus low conditional power",
    "NCT01326312": "FDA Clinical Hold",
    "NCT01481116": "Potential liver safety (hepatotoxicity) concerns",
    "NCT02532010": "Placed on clinical hold and never reopened",
    "NCT00496223": "Toxicities required dose reductions compromising effectiveness",
    "NCT03720574": "Halted to reassess safety risks and benefits",
    "NCT00953576": "Concerns regarding drug toxicity",
    "NCT04952753": "Urgent safety measures and halt of NIS793 treatment",
    "NCT03663166": "Increased SAE occurrence per PI",
    "NCT02711137": "Terminated due to safety issues",
    "NCT05360264": "Unexpected toxicity events during interim safety evaluation",
    "NCT00277316": "XL999 development halted due to cardiac toxicities in subjects",
    "NCT00393380": "Stopped because of toxicity concerns",
    "NCT01342016": "Safety concern with active control drug",
    "NCT00366795": "Terminated due to safety reason",
    "NCT01029925": "DSMB closed study due to higher-than-expected risk/safety concerns",
    "NCT07087938": "FDA clinical hold due to insufficient preclinical safety/evidence",
    "NCT02910882": "Explicit safety termination",
    "NCT01576887": "IDMC recommendation for safety concerns",
    "NCT00290017": "FDA Clinical Hold (May 2007)",
    "NCT02361346": "Terminated following regulatory clinical hold",
    "NCT00417248": "Safety concerns of regimen plus negative ESCAPE sorafenib readout",
    "NCT03713944": "Extreme toxicity and thromboembolic events",
    "NCT03285646": "Urgent safety measure halted dosing and enrollment",
    "NCT00955240": "Unacceptable toxicity profile",
    "NCT00383487": "Extreme toxicity",
    "NCT00667277": "Unacceptable toxicity profile with lack of response activity",
    "NCT00788606": "Toxicity and efficacy data from parallel trial",
    "NCT00243204": "FDA Clinical Hold (May 2007)",
    "NCT05223647": "Data Safety Committee advised termination due to excess experimental toxicity",
    "NCT00277290": "Cardiac toxicities observed in study subjects",
    "NCT00199342": "Two dose-limiting toxicities at initial 60 mg/m2 dose level",
    "NCT01792635": "Premature termination due to safety concern",
    # Reclassified from Stratum 3 upon manual clinical inspection:
    "NCT01312935": "Cardiovascular safety signal: observations of blood pressure reductions (hypotension) in treated subjects",
}

ADMIN_95: dict[str, str] = {
    # Reclassified false-positive safety/futility cases from Strata 1 & 2:
    "NCT03288545": "Sponsor closed study after meeting enrollment and primary objectives for regulatory filing; explicitly not safety or futility",
    "NCT00401401": "MTD established and 16-month safety follow-up completed; administrative closure of extended 3-year follow-up",
    "NCT04623216": "Sponsor strategic decision; explicitly not based on any safety findings or safety concerns",
    "NCT03945227": "Enrollment difficulty due to competing trial and Korean Ministry of Food and Drug Safety (MFDS) approval of durvalumab",
    "NCT03694197": "Portfolio redundancy: long-term safety objective of Praluent already adequately evaluated in other studies",
    "NCT05589597": "Sponsor portfolio decision; explicitly unrelated to safety (EO4010 safety profile benign)",
    "NCT03765541": "Under-accrual closure after reaching protocol event-maturity milestone (33% OS deaths) enabling primary analysis",
    # Stratum 3 administrative, commercial, operational, and accrual stops:
    "NCT00330200": "Poor patient enrollment",
    "NCT01420666": "Withdrawn for administrative reason",
    "NCT02560038": "Slow enrollment and resource reallocation",
    "NCT01315873": "PI departed institution prior to completion",
    "NCT01264081": "Unable to reach stated patient accrual",
    "NCT00949988": "Only 3 subjects enrolled and all 3 withdrew consent",
    "NCT05329168": "Business decision prior to first-patient-first-visit; explicitly not safety",
    "NCT03355482": "Poor recruitment and insufficient funding",
    "NCT04995094": "Enrollment challenges; study operationally infeasible",
    "NCT05177770": "Terminated for strategic reasons; explicitly not due to safety concerns",
    "NCT02413827": "Feasibility concerns due to evolving standard of care",
    "NCT01570192": "NIAID terminated due to low subject enrollment",
    "NCT04965688": "Insufficient enrollment to meet study goals",
    "NCT00239720": "Operational impact of hOKT3gamma1 (Ala-Ala) manufacturer change",
    "NCT02812693": "Poor patient accrual",
    "NCT04905407": "Sponsor administrative decision",
    "NCT00920556": "Administrative closure after 24 enrolled subjects provided sufficient decision data",
    "NCT05470985": "Terminated due to corporate business reasons",
    "NCT02716688": "Withdrawn with zero patients enrolled",
    "NCT06608940": "Funding unavailable",
    "NCT00690885": "Insufficient patient accrual",
    "NCT00983346": "Required laboratory staff left institution; accrual closed",
    "NCT00885963": "Financial resource limitation",
    "NCT04802031": "Low patient accrual",
    "NCT04134208": "Slow patient accrual",
    "NCT02508636": "Low patient accrual",
    "NCT02978443": "Lack of patient enrollment",
    "NCT00109850": "Closed due to poor accrual",
    "NCT07090174": "Withdrawn prior to enrollment due to lack of study drug supply",
    "NCT03241927": "Difficult patient enrollment",
    "NCT00846391": "Inability to recruit patients",
    "NCT02980705": "Sponsor administrative decision",
    "NCT00356525": "Stopped early due to low enrollment",
    "NCT00184028": "Insufficient patient accrual",
    "NCT04249622": "Slow patient accrual",
    "NCT04678323": "Insufficient funding",
    "NCT01897454": "Slower-than-anticipated patient accrual",
    "NCT05245760": "Funding unavailable",
    "NCT05423197": "Operational/logistical issues",
    "NCT02153905": "Slow and insufficient accrual",
    "NCT00003532": "Slow patient accrual",
    "NCT00349752": "Recruitment halted due to slow enrollment; no safety concerns",
    "NCT00304200": "Inadequate funding; sponsor companies requested closure",
    "NCT00356200": "Administrative closure after enrollment criteria met",
    "NCT04639024": "Corporate resource reallocation to Phase III program",
    "NCT02127073": "Unable to secure collaborator or funding after institutional transfer",
    "NCT02405910": "PI departed institution",
    "NCT03287271": "Loss of funding",
    "NCT03496298": "Sponsor decision to cancel trial; explicitly unrelated to safety",
    "NCT00372996": "Closed due to business reasons; explicitly not prompted by safety or efficacy concerns",
    "NCT00215631": "Slow recruitment",
    "NCT02928029": "Standard-of-care shift and slow recruitment",
    "NCT04975516": "Low/slow patient accrual",
    "NCT00544076": "Accrual too slow",
    "NCT02484209": "Overly strict inclusion criteria impeded primary-care recruitment",
    "NCT03295552": "Slow enrollment",
    "NCT00809341": "Low accrual",
    "NCT04892446": "Sponsor corporate decision to discontinue program development",
    "NCT00515112": "Poor patient accrual",
    # Reclassified from Stratum 4 upon manual clinical inspection:
    "NCT04522375": "Zero eligible participants identified 20 months after site initiation",
    "NCT00346749": "Difficulties identifying participating sites and willing subjects",
    "NCT00664105": "Experimental treatment became standard of care",
    "NCT05647512": "Withdrawn according to sponsor corporate strategy",
    "NCT01060345": "Study drug supply no longer available",
    "NCT00418821": "Administrative closure after FDA acknowledged fulfillment of Post-Marketing Commitment",
    "NCT01367275": "Sponsor administrative study closure",
    "NCT00447291": "Manufacturer discontinued fish-oil supplement and study food supply",
    "NCT01569074": "Corporate licensing return: AstraZeneca returned fostamatinib rights to Rigel",
    "NCT00386399": "All screened subjects tested BRCA2-negative; zero eligible patients to initiate trial",
    "NCT03360448": "Corporate clinical development strategy re-evaluation",
    "NCT03777124": "Sponsor R&D portfolio strategy adjustment",
    "NCT00426127": "Inadequate number of eligible patients for accrual",
    "NCT03038477": "Withdrawn/cancelled by principal investigator prior to enrollment",
    "NCT02701517": "Stopped prior to target enrollment due to lack of financial/institutional support",
    "NCT04732286": "Protocol-specified administrative closure after 24-month follow-up window",
    "NCT00493428": "Withdrawn prior to enrollment for protocol modification",
    "NCT01067326": "Sponsor-wide corporate termination of aliskiren trial portfolio",
    "NCT06191250": "Terminated prior to activation due to shift in standard of care",
    "NCT03858972": "Sponsor corporate discontinuation of tesetaxel development",
    "NCT04171284": "Corporate strategy adjustment; explicitly independent of safety and efficacy",
    "NCT03365661": "Withdrawn per sponsor decision",
    "NCT04453371": "Local ethics/regulatory committee denied approval request prior to enrollment",
    "NCT02713828": "Company discontinued commercial development for this indication",
    "NCT05564338": "Company-level development strategy change; explicitly not due to safety or health authority request",
    "NCT03228186": "Pharmaceutical company discontinued study drug supply/development",
    "NCT00242372": "Sponsor terminated clinical development program",
    "NCT00539799": "Local pharmacy operationally unwilling to comply with study protocol",
    "NCT04216472": "Insufficient participation rate (<75% participation)",
    "NCT01966575": "Withdrawn prior to start due to external literature favoring letrozole over clomiphene",
}

OTHER_UNREPORTED_26: dict[str, str] = {
    "NCT01339910": "DSMB recommended accrual termination without disclosing futility vs safety vs accrual cause",
    "NCT00243035": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT01097083": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT00264368": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT02773797": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT03906071": "Ambiguous registry text ('Inability to meet protocol objectives' without specifying accrual vs efficacy)",
    "NCT00325455": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT00347165": "Uninformative registry text ('undefined')",
    "NCT02421354": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT00051688": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT00073450": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT01484860": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT01866423": "Ambiguous registry text ('Not progressing toward scientific goals' without specifying accrual vs efficacy)",
    "NCT01501370": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT00695604": "Uninformative registry text ('Study was not able to be completed, no results analyzed')",
    "NCT00042302": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT00293436": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT03576885": "Unspecified DSMB recommendation at interim analysis without disclosed direction or cause",
    "NCT02458963": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT00188045": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT01884675": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT03343977": "Ambiguous registry text ('Study withdrawn for technical reasons')",
    "NCT02385669": "Placeholder registry text ('why study stopped')",
    "NCT00498901": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT00273507": "Unreported stop reason (blank why_stopped field on ClinicalTrials.gov)",
    "NCT01135329": "Unspecified stopping rule met without disclosing efficacy, futility, or toxicity direction",
}


def main() -> None:
    sample_p = ROOT / "outputs" / "stage2_tables" / "Table_S2a_Stoppage_Audit_240_Sample.csv"
    df = pd.read_csv(sample_p)
    assert len(df) == 240, f"Expected 240 trials, got {len(df)}"

    manual_map: dict[str, tuple[str, str]] = {}
    for cls_name, d in [
        ("Scientific_Efficacy_Futility", EFFICACY_59),
        ("Scientific_Safety_Toxicity", SAFETY_54),
        ("Administrative_Commercial_Enrollment", ADMIN_95),
        ("Other_or_Unreported", OTHER_UNREPORTED_26),
    ]:
        for nct, note in d.items():
            if nct in manual_map:
                raise ValueError(f"Duplicate NCT {nct} in manual adjudication dictionaries!")
            manual_map[nct] = (cls_name, note)

    missing = set(df["nct_id"]) - set(manual_map)
    extra = set(manual_map) - set(df["nct_id"])
    if missing or extra:
        raise ValueError(f"NCT mismatch! missing={missing}, extra={extra}")

    ledger_rows = []
    for row in df.itertuples(index=False):
        lbl, rationale = manual_map[row.nct_id]
        ledger_rows.append(
            {
                "nct_id": row.nct_id,
                "disease_axis": row.disease_axis,
                "overall_status": row.overall_status,
                "why_stopped": "" if pd.isna(row.why_stopped) else str(row.why_stopped),
                "manual_clinical_label": lbl,
                "manual_adjudication_rationale": rationale,
            }
        )

    out_df = pd.DataFrame(ledger_rows)
    out_p = ROOT / "data_vault" / "raw_clinical_trials" / "stoppage_240_manual_clinical_adjudication_ledger.csv"
    out_df.to_csv(out_p, index=False, encoding="utf-8")
    print(f"[Manual Adjudication Ledger Saved] {len(out_df)} trials -> {out_p}")
    print(out_df["manual_clinical_label"].value_counts().to_string())


if __name__ == "__main__":
    main()
