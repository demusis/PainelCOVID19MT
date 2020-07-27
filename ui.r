library(ggiraph)
library(shiny)
library(shinythemes)
library(stringi)
library(waiter)

# setwd("D:/OneDrive/Shiny")

# Carrega as informações por município.
aglomerados <-
  read.csv('aglomerados.csv',
           header = TRUE,
           stringsAsFactors = FALSE)
# Apaga eventuais celulas vazias.
aglomerados <- aglomerados[!is.na(aglomerados$aglomerado),]

l_aglomerados <- levels(factor(aglomerados$aglomerado))
l_aglomerados <- c('Mato Grosso', l_aglomerados)

# Define UI. 
shinyUI(
  fluidPage(theme = shinytheme('sandstone'),
            use_waiter(), # dependencies
            waiter_show_on_load(tagList(spin_heartbeat(),
                                        br(), br(), "Carregando..."
                                       )
            ),         
            
  tittle =  'Gabinete Militar: Análise das séries temporais da COVID-19 do Estado de Mato Grosso',
  
  # Título da aplicação
  
  titlePanel(title=div(img(src="brasao.png",
                           height = 50,
                           width = 50,
                           style = "margin:15px 15px"), 
                       "Gabinete Militar"),
             windowTitle = 'Gabinete Militar: Análise das séries temporais da COVID-19 do Estado de Mato Grosso'
             ),
  HTML(paste(
    p(strong(stri_dup(intToUtf8(160), 27),  
             "Análise das séries temporais da COVID-19 do Estado de Mato Grosso.")
  ))),
  
  
  # Sidebar com controles e previsões.
  sidebarPanel(
    width = 3, 
    fileInput('arq_aglomerados', 'Arquivos de dados:',
              multiple = TRUE,
              buttonLabel = 'Carregar',
              placeholder = 'Nenhum arquivo selecionado', 
              accept = c("csv", ".csv")),
    
    # linha horizontal ----
    tags$hr(),
    
    selectInput("aglomerado", "Aglomerado:",
                setNames(l_aglomerados, l_aglomerados)),
    selectInput("varepi", "Variável:",
                list('Infectados'= 'Infectados', 'Recuperados'='Recuperados', 'Óbitos'='Obitos', 'Total'='Total')),
    numericInput("ahead", "Número de registros a estimar:", 5),
    
    submitButton("Atualiza", icon("refresh"))
  ),
  

  
  # Apresenta a legenda e os gráficos de previsões.
  mainPanel(
    h3(textOutput("caption")),
    tabsetPanel(
      tabPanel("Média móvel", plotOutput("mmPlot"),
                              dataTableOutput('mmTabela'),
                              downloadButton('mmArquivo')),
      tabPanel("Bagged", plotOutput("baggedForecastPlot"),
                         dataTableOutput('baggedTabela'),
                         downloadButton('baggedArquivo')), 
      tabPanel("ARIMA", plotOutput("arimaForecastPlot"),
                        plotOutput("decomposicaoPlot"),
                        dataTableOutput('arimaTabela'),
                        downloadButton('arimaArquivo'),),
      tabPanel("Número básico de reprodução", textOutput("rtTexto"),
                                              plotOutput("mm7rtPlot"),
                                              plotOutput("rtPlot")),
      tabPanel("Velocidade de avanço", plotOutput("vaPlot")),
      tabPanel("Wavelet", plotOutput("waveletPlot")),
      tabPanel("Mapa de risco (lógica nebulosa)", girafeOutput("fuzzyPlot",  width = '800px', 
                                                                           height = '800px'),
                                                  dataTableOutput('fuzzyTabela'),
                                                  downloadButton('dwArquivo'),),
      tabPanel("Picos por ERS/Estado", dataTableOutput("sirTabela"))
    )
  )
))
