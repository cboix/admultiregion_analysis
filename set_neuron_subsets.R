#!/usr/bin/R
# ----------------------------------------
# Excitatory neuron to region annotations:
# Updated 03/10/2022 
# ----------------------------------------

ec.exc.neurons = c("Exc TOX3 POSTN", "Exc DLC1 SNTG2", "Exc AGBL1 GPC5", 
                   "Exc RELN GPC5", "Exc RELN COL5A2", "Exc TOX3 TTC6", 
                   "Exc TOX3 INO80D", "Exc SOX11 NCKAP5", "Exc COL25A1 SEMA3D")

hc.exc.neurons = c("DG granule cells", "CA1 pyramidal cells", 
                   "CA2, CA3 pyramidal cells", "Exc COBLL1 UST",
                   "Exc TRPC6 ANO2", "Exc GRIK1 CTXND1 (Subiculum)", 
                   "Exc ZNF385D COL24A1")

th.exc.neurons = c("Exc NXPH1 RNF220", "Exc VAT1L ERBB4",
                   "Exc SV2C LINC02137")

ctx.exc.neurons = c("Exc L2-3 CBLN2 LINC02306", "Exc L3-4 RORB CUX2",
                    "Exc L3-5 RORB PLCH1", "Exc L4-5 RORB GABRG1", 
                    "Exc L4-5 RORB IL1RAPL2", "Exc L5/6 IT Car3",
                    "Exc L5/6 NP", "Exc L5-6 RORB LINC02196",
                    "Exc L5 ET", "Exc L6b",
                    "Exc L6 CT", "Exc L6 THEMIS NFIA")

exc.set.regions = list('HCneurons'='HC',
                       'ECneurons'='EC',
                       'THneurons'='TH',
                       'CTXneurons'=c('AG','MT','PFC'))
exc.sets =  list('HCneurons'=hc.exc.neurons,
                 'ECneurons'=ec.exc.neurons,
                 'THneurons'=th.exc.neurons,
                 'CTXneurons'=ctx.exc.neurons)

vuln.neurons = c("Exc AGBL1 GPC5", "Exc RELN GPC5", 
                 "Exc RELN COL5A2", "Exc TOX3 TTC6", 
                 "CA1 pyramidal cells")

