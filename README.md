- 👋 Hi, I’m Charlotte Jeschina, Medical Student at the Unniversity of Lübeck and PhD-Student in the Network Transfer Study @Network-Transfer-HL

These are the scripts for my fMRI Analysis - Credits are at the tops of the scripts

This is the order in which I executed them:
# 1st-Level-Analysis:
- NetzTran_Timing_Files_MEMORY_Immediate_Recall_v5.m
- NetzTran_Timing_Files_MEMORY_Delayed_Recall_v2.m
- my_fMRI_FirstLevelAnalysis_concatenated_ImmediateRecall.m
- my_fMRI_FirstLevelAnalysis_concatenated_Delayed_Recall.m
- contrasts_1st_level_analysis_IR_DR_v2.m
Diese Skripte wurden verwendet, um:
- Ereigniszeiten zu definieren (Timing-Files),
- das statistische Modell zu erstellen (First-Level-Analyse),
- Kontraste zu definieren, z. B. IR > DR etc.
# 2nd-Level-Analysis:
- onesamplettest_fMRI_2ndLevel_FPA_DR.m
- onesamplettest_fMRI_2ndLevel_FPA_IR.m
- onesamplettest_fMRI_2ndLevel_NSWP_DR.m
- onesamplettest_fMRI_2ndLevel_NSWP_IR.m
- Anova_IQ_cov2_fMRI_2ndLevel.m
Diese Skripte dienten der Gruppenvergleichsanalyse der Aktivierungen sowie der Kovariatenkontrolle
