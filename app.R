# install.packages(c("shiny","DT"))
library(shiny)
library(DT)

nonOH2_comp <- c("CDCA","LCA","UDCA","GCDCA","GLCA","GUDCA","TCDCA","TLCA","TUDCA")
PriBA_comp  <- c("CA","CDCA","GCA","GCDCA","TCA","TCDCA")

taxa_up   <- c("GDCA", "GLCA", "TDCA", "GLCA_3S", "isoLCA", "TDCA_p")  # BAM2
taxa_down <- c("nonOH2_p", "PriBA_p")                                   # BAM1

input_ba <- unique(c(nonOH2_comp, PriBA_comp, "GDCA","GLCA","TDCA","GLCA_3S","isoLCA"))

# -------------------------- BDI calculator--------------------------
calc_BDI <- function(feature_df){
  mean_log_up   <- mean(log(as.numeric(feature_df[, taxa_up])), na.rm = TRUE)
  mean_log_down <- mean(log(as.numeric(feature_df[, taxa_down])), na.rm = TRUE)
  bdi <- mean_log_up - mean_log_down
  return(bdi)
}

# -------------------------- UI --------------------------
ui <- fluidPage(
  titlePanel("BDI Calculator | serum Bile Acid Dysmetabolism Index_FAH-SYSU"),
  hr(),
  # 横向输入区域，一行多个输入框，宽度缩短
  fluidRow(
    lapply(c("Total BA", input_ba), function(x){
      column(width = 2,
             numericInput(inputId = x, label = x, value = 0, min = 0, width = "100px")
      )
    })
  ),
  br(),
  actionButton("calc_btn", "Calculate BDI", class = "btn-primary"),
  helpText("Unit: μmol/L. Zero values will be substituted with 1e-6 to avoid log(0) error."),
  hr(),
  
  # 结果区域纵向依次排布
  h4("Calculation Log"),
  verbatimTextOutput("log_txt"),
  br(),
  
  h4("Result Table"),
  DTOutput("result_table"),
  br(),
  
  h4("Formula Notes"),
  p("nonOH2 = CDCA + LCA + UDCA + GCDCA + GLCA + GUDCA + TCDCA + TLCA + TUDCA"),
  p("PriBA = CA + CDCA + GCA + GCDCA + TCA + TCDCA"),
  p("nonOH2_p = nonOH2 / Total BA; PriBA_p = PriBA / Total BA; TDCA_p = TDCA / Total BA"),
  p("BDI = mean(log(GDCA, GLCA, TDCA, GLCA_3S, isoLCA, TDCA_p)) − mean(log(nonOH2_p, PriBA_p))"),
  br(),
  
  h4("Abbreviation List"),
  p("BA, bile acids."),
  p("CA, Cholic acid; CDCA, Chenodeoxycholic acid; GCA, Glycocholic acid; GCDCA, Glycochenodeoxycholic acid; GDCA, glycodeoxycholic acid; GLCA, Glycolithocholic acid; GUDCA, Glycoursodeoxycholic acid; LCA, Lithocholic acid; TCA, Taurocholic acid; TCDCA, Taurochenodesoxycholic acid; TDCA, Taurodeoxycholic acid; TLCA, Taurolithocholic acid; TUDCA, Tauroursodeoxycholic acid; UDCA, Ursodeoxycholic acid; GLCA_3S, Glycolithocholic acid 3-sulfate; isoLCA, Isolithocholic acid."),
  p("nonOH2, non-12α-Hydroxylated BAs and their conjugated forms. PriBA, sum of primary BAs."),
  br(),
  hr(),
  
  # Disclaimer 放在页面最底部
  p(strong("Disclaimer: "), "For research use only, not clinical use."),
  p("The First Affiliated Hospital, Sun Yat-Sen University (FAH-SYSU), 58 Zhongshan Er Road, Guangzhou 510080, China.")
)

# -------------------------- Server --------------------------
server <- function(input, output) {
  
  observeEvent(input$calc_btn, {
    ba_raw <- sapply(input_ba, function(x) as.numeric(input[[x]]))
    TotalBA <- as.numeric(input$TotalBA)
    
    ba_raw[ba_raw <= 0] <- 1e-6
    nonOH2 <- sum(ba_raw[nonOH2_comp])
    PriBA  <- sum(ba_raw[PriBA_comp])
    
    if(TotalBA <= 0){
      output$log_txt <- renderPrint({"Error: Total BA must be greater than zero!"})
      return()
    }
    
    nonOH2_p <- nonOH2 / TotalBA
    PriBA_p  <- PriBA / TotalBA
    TDCA_p   <- ba_raw["TDCA"] / TotalBA
    feature_df <- data.frame(
      GDCA     = ba_raw["GDCA"],
      GLCA     = ba_raw["GLCA"],
      TDCA     = ba_raw["TDCA"],
      GLCA_3S  = ba_raw["GLCA_3S"],
      isoLCA   = ba_raw["isoLCA"],
      TDCA_p   = TDCA_p,
      nonOH2_p = nonOH2_p,
      PriBA_p  = PriBA_p
    )
    
    bdi_score <- calc_BDI(feature_df)
    
    log_text <- paste0(
      "User-input Total BA = ", round(TotalBA,4), "\n",
      "nonOH2 = ", round(nonOH2,4), " | nonOH2_p = ", round(nonOH2_p,4), "\n",
      "PriBA = ", round(PriBA,4), " | PriBA_p = ", round(PriBA_p,4), "\n",
      "TDCA = ", round(ba_raw["TDCA"],4), " | TDCA_p = ", round(TDCA_p,4), "\n",
      "BDI score = ", round(bdi_score,4)
    )
    output$log_txt <- renderPrint(cat(log_text))
    
    res_df <- data.frame(
      Feature = colnames(feature_df),
      Value = round(as.numeric(feature_df),4)
    )
    res_df <- rbind(res_df, data.frame(Feature="BDI_score", Value=round(bdi_score,4)))
    output$result_table <- renderDT({
      datatable(res_df, rownames=F, options = list(dom='t'))
    })
  })
}

shinyApp(ui = ui, server = server)
