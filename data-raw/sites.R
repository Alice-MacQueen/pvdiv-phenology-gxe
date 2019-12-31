library(tidyverse)

sites <- tribble(~Site, ~Latitude, ~Longitude,
                 "PKLE", 30.383979, -97.729383,
                 "FRMI", 41.837133, -88.239667,
                 "OVTN", 32.3029, -94.9794,
                 "KING", 27.549863, -97.881012,
                 "BRKG", 44.3068, -96.6705,
                 "CLMB", 38.8969, -92.2178,
                 "KBSM", 42.419618, -85.371266,
                 "LINC", 41.1543, -96.4153,
                 "STIL", 35.991148, -97.046489,
                 "TMPL", 31.043383, -97.349498
)

save(sites, file = "data/sites.rda")
