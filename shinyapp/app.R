# make shinyapp of model correlograms for a set of phylogenetic models
# created 2021-02-02

library(shiny)

aalist = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V") # alpha by AA name

#aalist = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") # alpha by letter
#aalist = c("D", "P", "E", "N", "K", "R", "Q", "S", "G", "H", "T", "A", "C", "Y", "M", "V", "W", "L", "I", "F") # hpi

# data from https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html
# which itself derives from:
# Relationship of sidechain hydrophobicity and a‐helical propensity on the stability of the single‐stranded amphipathic a‐helix
# https://doi.org/10.1002/psc.310010507
hydrophobicity = c(41,-14,-28,-55,49,-10,-31,0,8,99,97,-23,74,100,-46,-5,13,96,63,76)
# sort and return index, which will reorder from hydrophilic (Asp) to most hydrophobic (Phe)
s_hyd = sort(hydrophobicity, index.return=TRUE)

modeldir = "~/git/graphphylo/mixture_model/"
#
wag_model = "~/git/graphphylo/mixture_model/WAG_model.hpi.txt"
# LG model of Le and Gascuel 2008 Mol Biol Evo
lg_model = "~/git/graphphylo/mixture_model/LG_model.hpi.txt"
# EX2 model of Le, Lartillot, Gascuel 2008 Phil. Trans. R. Soc. B
ex2_exp_model = "~/git/graphphylo/mixture_model/exp_EX2_model.hpi.txt"
ex2_bur_model = "~/git/graphphylo/mixture_model/bur_EX2_model.hpi.txt"
# UL2 model
ul2_m1_model = "~/git/graphphylo/mixture_model/ul2_m1_model.hpi.txt"
ul2_m2_model = "~/git/graphphylo/mixture_model/ul2_m2_model.hpi.txt"
# EX_EHO model of Le and Gascuel 2010 Systematic Biology
exp_ext_model = "~/git/graphphylo/mixture_model/exp_ext_model.hpi.txt"
bur_ext_model = "~/git/graphphylo/mixture_model/bur_ext_model.hpi.txt"
exp_hel_model = "~/git/graphphylo/mixture_model/exp_hel_model.hpi.txt"
bur_hel_model = "~/git/graphphylo/mixture_model/bur_hel_model.hpi.txt"
exp_oth_model = "~/git/graphphylo/mixture_model/exp_oth_model.hpi.txt"
bur_oth_model = "~/git/graphphylo/mixture_model/bur_oth_model.hpi.txt"


ui <- fluidPage(
    
    # App title ----
    titlePanel("Mixture models in phylogenetics"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            selectInput("phylomodel", h3("Phylogenetic Model"), 
                        choices = list("LG" = lg_model, "WAG" = wag_model, 
                                       "EX2 Exposed" = ex2_exp_model, "EX2 Buried" = ex2_bur_model,
                                       "UL2 component M1" = ul2_m1_model, "UL2 component M2" = ul2_m2_model,
                                       "EX_EHO exposed beta-sheet" = exp_ext_model, "EX_EHO buried beta-sheet" = bur_ext_model,
                                       "EX_EHO exposed alpha-helix" = exp_hel_model, "EX_EHO buried alpha-helix" = bur_hel_model,
                                       "EX_EHO exposed other" = exp_oth_model, "EX_EHO buried other" = bur_oth_model
                        ), 
                        selected = 1
            ),
            radioButtons("color_set", h3("Color set"),
                         choices = list("Blues (default)" = "blue", "Brown to Teal" = "teal", 
                                        "Purple to Green" = "purple")
            )
        ),
        # Main panel for displaying outputs ----
        mainPanel( 
          #  h3("Each point is a contig. Click-and-drag to display stats"),
            strong( textOutput("chosen_model") ),
            strong( textOutput("chosen_color") ),
          #  strong( paste("Using", basename(selectInput) ) ),
            plotOutput(outputId = "phylocorrplot",
                       height="600px",
                       click = "plot_click"
            ),
            tableOutput("info")
        )
    )
)

# Define server logic ----
server <- function(input, output) {
    color_palettes = list("blue" = c("#7fcdbb99", "#41b6c4aa", "#1d91c0aa", "#225ea8cc", "#081d58ee"),
                          "teal" = c("#bf812d", "#f6e8c3", "#c7eae5", "#35978f", "#003c30"),
                          "purple" = c("#762a83", "#c2a5cf", "#d9f0d3", "#a6dba0", "#1b7837")
    )
    
    output$chosen_model <- renderText({ 
        paste("Showing", basename(input$phylomodel) )
    })
    output$chosen_color <- renderText({ 
        paste("Coloring in", input$color_set )
    })
    output$chosen_color <- renderText({ 
        paste(" ", input$color_set )
    })
    
    output$phylocorrplot <- renderPlot({
        modeldat = read.table(input$phylomodel, sep=" ", fill=TRUE, col.names=aalist[s_hyd$ix])
        par(mar=c(4,1,1,4))
        plot(0,0,type="n",xlim=c(0,20),ylim=c(1,20),axes=FALSE, xlab="", ylab="")
        axis(1,at=1:20,labels=FALSE, cex.axis=1.4)
        mtext(aalist[s_hyd$ix], side=1, at=1:20, cex=1.4, line=1)
        axis(4,at=1:20,labels=FALSE, cex.axis=1.4)
        mtext(aalist[s_hyd$ix], side=4, at=1:20, cex=1.4, line=1)
        for (i in 1:19) {
            pointscaling = (6 + log(as.numeric(modeldat[i,1:i])))/2
            colorindex = floor(pointscaling)+1
            colorindex[colorindex<1]=1
            points(rep(i+1,i), 1:i, cex=4, pch=15, col= color_palettes[[input$color_set]][colorindex])
        }
        legend(0,20,legend=c("unfavored transition", "", "", "", "favored transition"), col=color_palettes[[input$color_set]], bty='n', pt.cex=3, cex=1.8, pch=15 )
    })

}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
