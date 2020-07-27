library(baytrends)
library(dplyr)
library(rdrop2)
library(EpiEstim)
library(forecast)
library(ggplot2)
library(ggiraph)
library(imputeTS)
library(lubridate)
library(openxlsx)
library(RColorBrewer)
library(shiny)
library(smooth)
library(sf)
library(waiter)
library(WaveletComp) # http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf


token <<- readRDS("droptoken.rds")

# setwd("D:/OneDrive/covid19mt")
carregaDados <- function() {
  # Carrega os MMC SIR por municipio.
  drop_download(dtoken = token, 'covid19mt/aux_sim.csv', overwrite = TRUE)
  sir <-
    read.csv("aux_sim.csv",
             header = TRUE,
             stringsAsFactors = FALSE)  
  
  # Carrega as informações por município.
  drop_download(dtoken = token, 'covid19mt/aglomerados.csv', overwrite = TRUE)
  aglomerados <-
    read.csv("aglomerados.csv",
             header = TRUE,
             stringsAsFactors = FALSE)
  
  # Apaga eventuais celulas vazias.
  aglomerados <- aglomerados[!is.na(aglomerados$aglomerado), ]
  
  # Converte para data
  aglomerados$data <-
    as.Date(parse_date_time(aglomerados[["data"]], '%d%m%y'), format = '%d%m%y')
  
  # Carrega as series temporais da epidemia.
  drop_download(dtoken = token, 'covid19mt/evolucao_ers.CSV', overwrite = TRUE)
  
  aux_dados <- read.csv("evolucao_ers.CSV",
                        stringsAsFactors = FALSE)
  aux_dados$Data <- dmy(aux_dados$Data)
  aux_dados$Dia_Juliano <- yday(aux_dados$Data)
  aux_dados$Total <-
    aux_dados$Infectados + aux_dados$Recuperados + aux_dados$Obitos
  
  # Apaga eventuais celulas vazias.
  # mt_dados <- aux_dados[!is.na(aux_dados[, 'Infectados']),]
  mt_dados <- aux_dados

  # Guarda data do último registro
  aux_dh <<- aux_dados[nrow(aux_dados), 'Data']
    
  # smt_aux_dados <- mt_dados[, c('Total', 'Municipio', 'Dia_Juliano')]
  
  # Agrega todos os dados para o Mato Grosso.
  mt_total_dados <-
    aggregate(cbind(Infectados, Recuperados, Obitos, Total) ~ Dia_Juliano,
              mt_dados,
              sum, na.action=na.pass)
  mt_total_dados$Municipio <- 'Mato Grosso'
  
  # Adiciona o Mato Grosso no rol.
  mt_aux_dados <<-
    rbind(mt_dados[, c('Municipio',
                       'Infectados',
                       'Recuperados',
                       'Obitos',
                       'Total',
                       'Dia_Juliano')],
          mt_total_dados)
  
  # Cria fatores
  mt_aux_dados$Municipio <<- factor(mt_aux_dados$Municipio)
}

ts_infectados <- function(aux_nivel, aux_varepi) {
  municipio_dados <-
    mt_aux_dados[mt_aux_dados$Municipio == aux_nivel, ]
  seq_dados <-
    data.frame(seq(
      min(municipio_dados$Dia_Juliano),
      max(municipio_dados$Dia_Juliano),
      by = 1
    ))
  colnames(seq_dados) <- c('Dia_Juliano')
  
  # Cria serie numerica fechada para verificar ausentes.
  todos_dados <-
    merge(x = seq_dados,
          y = municipio_dados,
          by = 'Dia_Juliano',
          all.x = TRUE)
  todos_dados$Municipio <- aux_nivel
  
  todos_dados$Infectados <- fillMissing(todos_dados$Infectados, span=2, max.fill = 90)
  todos_dados$Recuperados <- fillMissing(todos_dados$Recuperados, span=2, max.fill = 90)
  todos_dados$Obitos <- fillMissing(todos_dados$Obitos, span=2, max.fill = 90)
  todos_dados$Total <- fillMissing(todos_dados$Total, span=2, max.fill = 90)
  
  # Seleciona variaveis de analise.
  # dados <- todos_dados[, c('Dia_Juliano', 'Total')]

  # Cria objeto TS.
  ts_aux <-
    ts(todos_dados[, aux_varepi],
       start = todos_dados[1, 'Dia_Juliano'],
       frequency = 1)
  
  # Preenchimento de falhas mais rigoroso.
  # try(ts_aux <-
  #   na_kalman(ts_aux, model = "StructTS", smooth = TRUE))

  return(ts_aux)
}

# Velocidade de avanço
va <- function(or_st) {
  st <- data.frame(dj = integer(),
                   y = double(),
                   stringsAsFactors = FALSE)
  posicao <- end(or_st)[1] - start(or_st)[1] + 1
  for (n_posicao in 15:posicao) {
    valor_t <- or_st[n_posicao]
    valor_t_1 <- or_st[n_posicao - 7][1]
    valor_t_2 <- or_st[n_posicao - 14][1]
    aux_st <- data.frame(time(or_st)[n_posicao],
                         (valor_t - valor_t_1) / (valor_t_1 - valor_t_2) - 1)
    colnames(aux_st) <- c('dj', 'y')
    st <- rbind(st, aux_st)
  }
  aux_st <- ts(st$y,
               start = st[1, 'dj'],
               frequency = 1)
  aux_st[is.nan(aux_st)] <- NA
  aux_st[is.infinite(aux_st)] <- NA
  return(aux_st)
}

# Mapa de risco.
drop_download(dtoken = token, 'covid19mt/fuzzy_mt.csv', overwrite = TRUE)
df <- read.csv(file = 'fuzzy_mt.csv')
municipios <- st_read("ERS_MUNICIPIOS_MT.shp")
aux_df <- df[, c('Geocodigo', 'risco')]
aux_df$Geocodigo <- as.character(aux_df$Geocodigo)
municipios <- municipios %>% left_join(aux_df, by = "Geocodigo")
tabela_municipios <- data.frame(municipios)[,c(1, 2, 12, 17, 18)]

# SIR
total_SIR <- function() {
  drop_download(dtoken = token, 'covid19mt/aux_sim.csv', overwrite = TRUE)
  sir <-
    read.csv("aux_sim.csv",
             header = TRUE,
             stringsAsFactors = FALSE)
  
  
  # Agrega todos os casos por dia e o associa ao MT.
  mt_df_sim <-
    aggregate(sir[, c('time',
                      's.num',
                      'infectados_oficial',
                      'recuperados_oficial')], list(sir$time), sum)[-2]
  
  mt_df_sim$ERS <- 'Mato Grosso'
  colnames(mt_df_sim) <-
    c('Data',
      'Suscetiveis',
      'Infectados',
      'Recuperados/Fatalidades',
      'ERS')
  mt_df_sim$ERS <- factor(mt_df_sim$ERS)
  mt_df_sim$Data <- aux_dh + mt_df_sim$Data - 1
  
  # Agrega os casos por regiao (ERS) e dia.
  aux_df_sim <-
    aggregate(sir[, c('time',
                      's.num',
                      'infectados_oficial',
                      'recuperados_oficial')],
              list(sir$regiao,
                   sir$time),
              sum)[-3]
  colnames(aux_df_sim) <-
    c('ERS',
      'Data',
      'Suscetiveis',
      'Infectados',
      'Recuperados/Fatalidades')
  aux_df_sim$ERS <- factor(aux_df_sim$ERS)
  aux_df_sim$Data <- aux_dh + aux_df_sim$Data - 1
  
  # Localiza ponto de maximo no Mato Grosso e o armazena
  ers_aux_df_sim <<- mt_df_sim
  ers_aux_df_sim <<-
    ers_aux_df_sim[which.max(ers_aux_df_sim$Infectados), ]
  
  # Localiza os maximos de infectados nas ERS e os armazena.
  for (aux_ers in levels(aux_df_sim$ERS)) {
    aux_ers_df_sim <- aux_df_sim[aux_df_sim$ERS == aux_ers,]
    aux_ers_df_sim <-
      aux_ers_df_sim[which.max(aux_ers_df_sim$Infectados), ]
    ers_aux_df_sim <<- rbind(ers_aux_df_sim, aux_ers_df_sim)
  }
}


#
# carrega os dados
#
carregaDados()
total_SIR()

# Inicia o servidor.
shinyServer(function(input, output) {
  waiter_hide()
  w <- Waiter$new(id = c('rtPlot', 'mm7rtPlot',
                         'baggedForecastPlot',
                         'arimaForecastPlot',
                         'decomposicaoPlot',
                         'arimaTabela',
                         'vaPlot',
                         'waveletPlot',
                         'fuzzyPlot'),
                  html = spin_gauge(),
                  color = 'white')
  
  getDataset <- reactive({
    return(ts_infectados(input$aglomerado, input$varepi))
  })
  
  obtemTotal <- reactive({
    return(ts_infectados(input$aglomerado, 'Total'))
  })
  
  atualizaDados <- reactive({
    if (!is.null(input$arq_aglomerados)) {
      carregaDados()
      nt <- showNotification("Dados carregados.",
                             duration = 5,
                             closeButton = FALSE,
                             type= 'warning')
      # Inserir procedimentos para verificar os dados e gravá-los em um arquivo permanente.
      write.csv(file = 'mt_aux_dados.csv', mt_aux_dados)
    }
  })  
  
  output$caption <- renderText({
    atualizaDados()
    paste('Último registro: ', aux_dh)
  })
  
  output$mmPlot <- renderPlot({
    w$show()
    nt <- showNotification("Processando médias móveis (7 dias).",
                           duration = NA,
                           closeButton = FALSE,
                           type= 'message')
    dados <- tail(diff(getDataset()), 30)
    mm <- sma(
      dados,
      order = 7,
      h = 1,
      interval = 'none',
      silent = "none"
    )
    removeNotification(req(nt))
  })
  
  output$mmTabela <- renderDataTable({
    nt <- showNotification("Processando tabelas das médias móveis.",
                           duration = 3,
                           closeButton = FALSE,
                           type= 'message')
    dados <- tail(diff(getDataset()), 30)
    mm <<- sma(
      dados,
      order = 7,
      h = 1,
      interval = 'none',
      silent = "all"
    )
    aux_tab <- cbind(mm$y, mm$fitted)
    aux_tab <- data.frame(aux_tab) 
    colnames(aux_tab) <- c('Indivíduos novos', 'Média móvel')
    aux_tab
  })
  
  output$mmArquivo <- downloadHandler(filename = 'mm.xlsx', 
                                      content = function(arqu){
                                         write.xlsx(mm, arqu)
                                      }
  )
  
  
  output$rtPlot <- renderPlot({
    w$show()
    nt <- showNotification("Processando número básico de reprodução.",
                           duration = NA,
                           closeButton = FALSE,
                           type= 'message')
    dados_MT <- as.numeric(diff(obtemTotal()))
    erro_rt <- try({
      res <- estimate_R(dados_MT, method = "parametric_si",
                        config = make_config(list(mean_si = 3.96, std_si = 4.75)))
      r30 <- tail(res$R$'Median(R)', 30)
      mm_r30 <- sma(
        r30,
        order = 7,
        h = 1,
        interval = 'none',
        silent = "none"
      )
      
    })
    removeNotification(req(nt))
  })
  
  output$mm7rtPlot <- renderPlot({
    w$show()
    nt <-
      showNotification(
        "Processando média movel do número básico de reprodução.",
        duration = 1,
        closeButton = FALSE,
        type = 'message'
      )
    dados_MT <- as.numeric(diff(obtemTotal()))
    erro_rt <- try({
      res <- estimate_R(dados_MT, method = "parametric_si",
                        config = make_config(list(mean_si = 3.96, std_si = 4.75)))
      plot(res)
    })
  })
  
  output$rtTexto <- renderText({
    nt <- showNotification("Processando intervalo de confiança do número básico de reprodução.",
                           duration = 1,
                           closeButton = FALSE,
                           type= 'message')
    dados_MT <- as.numeric(diff(obtemTotal()))
    erro_rt <- try({
      res <- estimate_R(dados_MT, method = "parametric_si",
                        config = make_config(list(mean_si = 3.96, std_si = 4.75)))
      r_025 <- round(tail(res$R$'Quantile.0.025(R)', 1), 3)
      r_975 <- round(tail(res$R$'Quantile.0.975(R)', 1), 3)
      r_500 <- round(tail(res$R$'Median(R)', 1), 3)
      paste('Intervalo de confiança da mediana (95%): [',
            r_025,
            ',',
            r_500,
            ',',
            r_975,
            '].')
    })
  })
  
  output$decomposicaoPlot <- renderPlot({
    w$show()
    nt <- showNotification("Processando decomposição da série.",
                           duration = 1,
                           closeButton = FALSE,
                           type= 'message')
    try({
      ts_aux <-
        ts(getDataset(), frequency = 7)
      f <- decompose(ts_aux)
      plot(f)
    })
  })

  output$arimaTabela <- renderDataTable({
    w$show()
    nt <- showNotification("Processando previsões do ARIMA.",
                           duration = 3,
                           closeButton = FALSE,
                           type = 'message')
    try({
      l <- BoxCox.lambda(getDataset())
      fit <- auto.arima(getDataset(), lambda = l)
      f <- forecast(fit, h = input$ahead)
      df_f <- data.frame(f)
      df_f <- df_f[, c(1, 4, 5)]
      colnames(df_f) <-
        c('Previsao', 'Limite inferior (95%)', 'Limite superior (95%)')
      df_arima <<- round(df_f, 2)
      df_arima
    })
  })
  
  output$arimaArquivo <- downloadHandler(filename = 'arima_estimativas.xlsx', 
                                         content = function(arqu){
                                           write.xlsx(df_arima, arqu)
                                         }
  )

  output$arimaForecastPlot <- renderPlot({
    w$show()
    nt <- showNotification("Processando gráfico do ARIMA.",
                           duration = 3,
                           closeButton = FALSE,
                           type= 'message')
    try({
      l <- BoxCox.lambda(getDataset())
      fit <- auto.arima(getDataset(), lambda = l)
      f <- forecast(fit, h = input$ahead)
      plot(forecast(f))
    })
  })
  
  output$baggedForecastPlot <- renderPlot({
    w$show()
    
    nt <- showNotification("Processando gráfico do Bagged.",
                           duration = 3,
                           closeButton = FALSE,
                           type= 'message')
    try({
      l <- BoxCox.lambda(getDataset())
      fit <- baggedModel(getDataset(), lambda = l)
      plot(forecast(fit, h = input$ahead))
    })
  })
  
  output$baggedTabela <- renderDataTable({
    nt <- showNotification("Processando tabelas do Bagged.",
                           duration = 3,
                           closeButton = FALSE,
                           type= 'message')
    try({
      l <- BoxCox.lambda(getDataset())
      fit <- baggedModel(getDataset(), lambda = l)
      f <- forecast(fit, h = input$ahead)
      df_f <- data.frame(f)
      colnames(df_f) <-
        c('Previsao', 'Limite inferior (95%)', 'Limite superior (95%)')
      df_bagged <<- round(df_f, 2)
      df_bagged
    })
  })
  
  output$baggedArquivo <- downloadHandler(filename = 'bagged_estimativas.xlsx', 
                                          content = function(arqu){
                                            write.xlsx(df_bagged, arqu)
                                          }
  )
  
  output$vaPlot <- renderPlot({
    w$show()
    nt <- showNotification("Processando média móvel da velocidade de avanço.",
                           duration = 1,
                           closeButton = FALSE,
                           type= 'message')
    try({
      ts_in1s <- va(getDataset())
      modelo_in1s <- sma(
        tail(ts_in1s, 14),
        order = 7,
        h = 1,
        interval = 'none',
        silent = 'none'
      )
    })
  })
  
  output$waveletPlot <- renderPlot({
    w$show()
    nt <- showNotification("Processando Wavelet.",
                           duration = NA,
                           closeButton = FALSE,
                           type= 'message')
    df_ts <- data.frame(getDataset())
    colnames(df_ts) <- 'x'
    try({
      my.w <- analyze.wavelet(
        df_ts,
        "x",
        loess.span = 0,
        dt = 1,
        dj = 1 / 100,
        lowerPeriod = 7,
        upperPeriod = 91,
        make.pval = TRUE,
        n.sim = 50
      )
      
      wt.image(
        my.w,
        color.key = "quantile",
        n.levels = 250,
        legend.params = list(lab = "wavelet power levels", mar = 4.7)
      )
    })
    removeNotification(nt)
  }, height = 700, width = 970)
  
  output$fuzzyPlot <- renderGirafe({
    w$show()
    nt <- showNotification("Processando mapa do risco por lógica fuzzy.",
                           duration = 3,
                           closeButton = FALSE,
                           type= 'message')
    aux_gg <- ggplot(municipios) +
      geom_sf_interactive(aes(fill = risco, tooltip = Nome), show.legend = "point") +
      theme_minimal() + 
      scale_fill_gradientn(colours=rev(brewer.pal(11, "Spectral")), 
                           na.value = "#ffffff",
                           limits = c(0, 100))
    girafe(ggobj = aux_gg, 
           options = list(opts_selection(type = "single", only_shiny = FALSE)))
  })   
  
  output$fuzzyTabela <- renderDataTable({
    nt <- showNotification("Processando tabela risco por lógica fuzzy.",
                           duration = 1,
                           closeButton = FALSE,
                           type= 'message')
    tabela_municipios
  })
  
  output$dwArquivo <- downloadHandler(filename = 'fuzzy_risco.xlsx', 
                                      content = function(arqu){
                                         write.xlsx(tabela_municipios, arqu)
                                      }
  )
  
  output$sirTabela <- renderDataTable({
    nt <- showNotification("Processando SIRs por ERS e Estado.",
                           duration = 2,
                           closeButton = FALSE,
                           type= 'message')
    ers_aux_df_sim
  })
})








