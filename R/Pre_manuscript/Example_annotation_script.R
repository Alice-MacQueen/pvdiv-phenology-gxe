library(tibble)
library(switchgrassGWAS)
library(AnnotationDbi)
txdb <- loadDb(file = file.path("/", "home", "SharedUser", "PV",
                                "Pvirgatum_516_v5.1.gene.txdb.sqlite"))

df <- tibble(CHR = "Chr03K", start = 50800000, end = 50900000)

annos <- pvdiv_table_topsnps(df = df, type = "table", txdb = txdb)
# annos should be a dataframe containing one gene, Pavir.3KG393200
