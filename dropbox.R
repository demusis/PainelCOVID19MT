library(rdrop2)

setwd("D:/OneDrive/covid19mt")

# Cria token
# token <- drop_auth()
# saveRDS(token, "droptoken.rds")

# Carrega token.
token <- readRDS("droptoken.rds")

# Then pass the token to each drop_ function
drop_acc(dtoken = token)

drop_download('covid19mt/Teste.txt', overwrite = TRUE)
drop_upload(dtoken = token, 'Teste.txt', path = 'covid19mt')

