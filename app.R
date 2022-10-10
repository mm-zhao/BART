#set working directory to directory containing BART repository
#install the following packages:
library(shiny)
options(java.parameters = "-Xmx4g")
library(bartMachine) # v1.2.6
library(shinydashboard)
library(shinydashboardPlus)
library(caret)
library(CORElearn)

#### FUNCTION
get_probability = function(df){
  df$msin = factor(df$msin, levels = c(0,1))
  df$msin = as.numeric(levels(df$msin))[df$msin]
  df$brafmt = factor(df$brafmt, levels = c(0,1))
  df$brafmt = as.numeric(levels(df$brafmt))[df$brafmt]
  df$female = factor(df$female, levels = c(0,1))
  df$female = as.numeric(levels(df$female))[df$female]
  df$tstage = factor(df$tstage, levels = c(1,2,3,4))
  df$tstage = as.numeric(levels(df$tstage))[df$tstage]
  df$sitea = as.numeric(df$sitea)
  df$necrosisp = ifelse(df$necrosisp == 0, 0, 
                        ifelse(df$necrosisp < 10, 1,
                               ifelse(df$necrosisp < 20, 2,
                                      ifelse(df$necrosisp < 30, 3, 
                                             ifelse(df$necrosisp < 40, 4,
                                                    5 )))))
  
  bart = readRDS("bart_test.rds")
  preproc_cont = readRDS("preproc_cont.rds")
  preproc_immune = readRDS("preproc_immune.rds")
  pre = c("ageyr", "intNoNLN", "necrosisp", "tstage", "sitea", "presmoker", "intNoPLN")
  immune = "density_tumor_cd3pcd4pcd45rop"
  df[pre] = predict(preproc_cont, df[pre])
  df[immune] = predict(preproc_immune, df[immune])
  
  pred_prob = predict(bart, df)

  quant = ifelse(pred_prob < 0.753, "High risk",
                 ifelse(pred_prob < 0.883, "Intermediate risk",
                        ifelse(pred_prob <= 1, "Low risk")))

  return(list(pred_prob = round(pred_prob,3), quant = quant))
}

#### UI
ui <- dashboardPage(
  dashboardHeaderPlus(title = "BART Risk Predictor"),
  dashboardSidebar(collapsed = T, width = 180,
                   sidebarMenu( style = 'font-size: 120%',
                     menuItem("Model", tabName = "home", icon = icon("connectdevelop")),
                     menuItem("Variable Definitions", tabName = "vardef"),
                     menuItem("About", tabName = "about"))),
  dashboardBody(
    tabItems(
    tabItem(tabName = "home",
  # Sidebar panel for inputs ----
    column(3, style = 'font-size: 120%',
           
         selectInput("stage", "TNM Stage:", #stage
                     c("Stage II" = "0",
                       "Stage III" = "1")),
         
         numericInput("ageyr", "Age (years):", #age
              value = ""),
         
         selectInput("sitea", "Tumor site:", #tumor depth
                     c("", "Cecum" = "186.25",
                       "Ascending colon" = "171.35",
                       "Hepatic flexure" = "159.8",
                       "Transverse colon" = "130.65",
                       "Splenic flexure" = "101.5",
                       "Descending colon" = "85",
                       "Sigmoid" = "44",
                       "Rectosigmoid" = "19.5",
                       "Rectal" = "9.75"))), 
  
    column(3,style = 'font-size: 120%',

    numericInput("intNoPLN", "No. of positive nodes:", #number of pln
                value = ""),
    
    selectInput("female", "Gender:", #gender
                c("",
                  "Male" = "0",
                  "Female" = "1")),
  
    numericInput("necrosisp", "Necrosis percentage:", #necrosis percentage
               value = "")),
  
    column(3, style = 'font-size: 120%',
           
    numericInput("intNoNLN", "No. of negative nodes:", #number of nln
                        value = ""),
           
     selectInput("msin", "MSI status:", #msi status
                 c("",
                   "MSI stable/MSI low" = "0",
                   "MSI high" = "1")),
    
    selectInput("brafmt", HTML("<i>BRAF</i> mutation:"), #BRAF mut
                c("", "Absent" = "0",
                  "Present" = "1"))),

  column(3, style = 'font-size: 120%',
         

    selectInput("tstage", "Tumor depth:", #tumor depth
                     c("", "T1" = "1",
                       "T2" = "2",
                       "T3" = "3",
                       "T4" = "4")),
  
    numericInput("presmoker", "Smoking history (pack-years):", #smoking hx
                 value = ""),
  
    numericInput("density_tumor_cd3pcd4pcd45rop", HTML("Tumor memory helper T cell density (cells/mm<sup>2</sup>):"), #tumor helper T density
                 value = "")),  
  
    column(11, 
           offset = 5,
              actionButton("do", "Run", style = "color: white; background-color: #59996a;  font-size:120%"),
              actionButton("reset", "Reset", style = "color: 
#848d95; background-color: #dce6df; font-size:120%")),
    column(12,
           HTML("<br> "),
           HTML("<br> ")),
  # Main panel for displaying outputs ----
  tags$head(tags$style(
    type="text/css",
    "#surv_q img {max-width: 100%; width: 100%; height: auto}",
    "#surv_curv img {max-width: 100%; width: 100%; height: auto}"
  )),
    fluidRow( style = 'font-size: 115%',
      column(6, offset = 0, style='padding:0px;', align = "center",
             gradientBox(
               title = "Results",
               gradientColor = "maroon",
               icon = "fa fa-bullseye",
               width = 12,
               footer = 
                 div(
                   htmlOutput("probability"),
               htmlOutput("quantile"))),
             gradientBox(
               title = "Risk groups by survival probability",
               icon = "fa fa-th",
               gradientColor = "light-blue",
               boxToolSize = "xs",
               width = 12,
               collapsed = T,
               footer =
                 imageOutput("surv_q", height = 200),
               "(Stage-specific)")
),
      column(6, offset = 0, style='padding:0px;', align = "center",
             gradientBox(
               title = "Survival curve",
               icon = "fa fa-signal",
               gradientColor = "teal", 
               boxToolSize = "sm",
               width = 12,
               collapsed = T,
               footer = 
                 imageOutput("surv_curv", height = 350),
               "(Stage-specific)")
)
      )
),

  tabItem(tabName = "vardef",
          h4("Definitions"),
          div(
          HTML("<b>Stage and tumor depth of invasion:</b>",
               "<br>Based on <a href=https://cancerstaging.org/references-tools/quickreferences/documents/colonmedium.pdf>AJCC Staging Criteria</a> for Colorectal Cancer.",
               "<br><br><b>Necrosis percentage:</b>",
               "<br>Estimated percentage of necrosis area noted within tumor on available whole tissue slides.",
               "<br><br><b>Number of negative LNs:</b>",
               "<br>The number of total lymph nodes discovered subtracted by the number of positive lymph nodes found."
))),
  tabItem(tabName = "about",
          h4("About"),
          div(
            HTML("BART Risk Predictor is a risk prediction model for <b>Stage II/Stage III Colorectal Cancer</b> based on an ensemble sum-of-trees learning model, Bayesian Additive Regression Trees, proposed by <a href = cite>Chipman et al.</a>"),
            HTML("<br>This model is trained on the <a href=https://www.nhs.us/>Nurses' Health Study (NHS)</a> and <a href = https://sites.sph.harvard.edu/hpfs/>Health Professionals Follow-Up Study (HPFS)</a> cohorts.",
                 "<br>For more information on the risk estimator, see <a href=http://tobepublished.com>here</a> for the orginial paper.")))
)
)
)


###### SERVER
server <- function(input, output, session) {
  observeEvent(input$reset, {
    updateNumericInput(session, "ageyr", value = "")
    updateNumericInput(session, "intNoNLN", value = "")
    updateNumericInput(session, "necrosisp", value = "")
    updateSelectInput(session, "msin", selected = "")
    updateSelectInput(session, "tstage", selected = "")
    updateNumericInput(session, "presmoker", value = "")
    updateNumericInput(session, "intNoPLN", value = "")
    updateSelectInput(session, "brafmt", selected = "")
    updateSelectInput(session, "female", selected = "")
    updateNumericInput(session, "density_tumor_cd3pcd4pcd45rop", value = "")
    updateSelectInput(session, "sitea", selected = "")
    updateSelectInput(session, "stage", selected = "")
    output$probability = renderText({})
    output$quantile = renderText({})
    output$surv_curv = renderText({})
    # output$surv_q = renderText({})
  })
  observeEvent(input$do, {
    df = data.frame(input$intNoNLN, input$tstage, input$ageyr, input$brafmt, input$msin, 
                    input$necrosisp, input$sitea, input$presmoker, input$intNoPLN, 
                    input$female, input$density_tumor_cd3pcd4pcd45rop)

    names(df) = c("intNoNLN", "tstage", "ageyr", "brafmt", "msin", "necrosisp", "sitea", "presmoker", "intNoPLN", "female", "density_tumor_cd3pcd4pcd45rop")
    pred = get_probability(df)
    stage = ifelse(input$stage == "0", "Stage II", "Stage III")

    
    if (pred$quant == "High risk"){
    output$probability = renderText({
      paste0("The estimated five-year survival probability is ",
             "<font color=\"#b5341d\"><b>",pred$pred_prob,"</b></font>",".")
    })
    output$quantile = renderText({
      paste0("This probability falls within risk group ", 
             "<font color=\"#b5341d\"><b>", pred$quant, "</b></font>",".")
    })
    } else if (pred$quant == "Intermediate risk"){
      output$probability = renderText({
        paste0("The estimated five-year survival probability is ",
               "<font color=\"#d6a62b\"><b>",pred$pred_prob,"</b></font>",".")
      })
      output$quantile = renderText({
        paste0("This probability falls within risk group ", 
               "<font color=\"#d6a62b\"><b>", pred$quant, "</b></font>", ".")
      })
    }else {
      output$probability = renderText({
        paste0("The estimated five-year survival probability is ",
               "<font color=\"#4e992b\"><b>",pred$pred_prob,"</b></font>",".")
      })
      output$quantile = renderText({
        paste0("This probability falls within risk group ", 
               "<font color=\"#4e992b\"><b>", pred$quant, "</b></font>", ".")
      })
    }
    output$surv_curv = renderImage({
        return(list(src='www/S2_and_3.png', align = "center", height = 250, width = 300))         
      }, deleteFile = FALSE)
    output$surv_q = renderImage({
        return(list(src='www/S2_and_3q.png', align = "center", width = 300))
      }, deleteFile = FALSE)
  })
}  

shinyApp(ui, server)
