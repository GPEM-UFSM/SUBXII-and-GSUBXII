#-------------------------------------------------
# Libraries
#-------------------------------------------------
library(shiny)
library(bslib)
library(ggplot2)
library(thematic)
library(shinycssloaders)
#-------------------------------------------------
#Enable thematic for automatic plot stylin
thematic_shiny(font = "auto")
#-------------------------------------------------
# MATH FUNCTIONS - Unit Sine Burr XII
#-------------------------------------------------
#Probability Distribution Function----
dSUBXII <- function(y, mu, sigma, log = FALSE) {
  if (any(mu <= 0) || any(sigma <= 0)) stop("mu and sigma must be positive")
  if (any(y <= 0 | y >= 1)) return(rep(NA, length(y)))
  fy <- -(pi * mu * sigma *
            (((-log(y))^sigma + 1)^mu + 2) *
            ((-log(y))^sigma + 1)^(-2 * mu - 1) *
            (-log(y))^sigma *
            cos((pi * (1 / ((-log(y))^sigma + 1)^mu + 1)) /
                  (4 * ((-log(y))^sigma + 1)^mu))
  ) / (4 * y * log(y))
  if (log) log(fy) else fy
}
#-------------------------------------------------
#Cumulative Distribution Function----
pSUBXII <- function(q, mu, sigma) {
  if (any(mu <= 0) || any(sigma <= 0)) stop("mu and sigma must be positive")
  if (any(q <= 0 | q >= 1)) return(rep(NA, length(q)))
  G <- (1 + (-log(q))^sigma)^(-mu)
  sin((pi / 4) * G * (G + 1))
}
#-------------------------------------------------
#Hazard Function----
hSUBXII <- function(y, mu, sigma) {
  f <- dSUBXII(y, mu, sigma)
  F <- pSUBXII(y, mu, sigma)
  f / (1 - F)
}
#-------------------------------------------------
#Survival Function----
sSUBXII <- function(y, mu, sigma) {
  1 - pSUBXII(y, mu, sigma)
}
#-------------------------------------------------
# UI
#-------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  withMathJax(),
  titlePanel("Unit Sine Burr XII"),
  hr(),
  tabsetPanel(
    tabPanel(title = "SUBXII Distribution",
             sidebarLayout(
               sidebarPanel(
                 h4("Parameters"),
                 sliderInput(inputId = "mu",
                             label=(helpText(c("$$\\mu$$"))),
                             min = 0.01, max = 20, value = 1.8, step = .2,
                             animate = animationOptions(interval = 2000, loop = T)),
                 sliderInput(inputId = "sigma",
                             label=(helpText(c("$$\\sigma$$"))),
                             min = 0.01, max = 20, value = 0.8, step = .2,
                             animate = animationOptions(interval = 2000, loop = T))),
               mainPanel(
                 fluidRow(
                   column(6, plotOutput("p1")), column(6, plotOutput("p2")),
                   column(6, plotOutput("p3")), column(6, plotOutput("p4"))
                 )
               )
             )
    )
  )
)
#-------------------------------------------------
# SERVER
#-------------------------------------------------
server <- function(input, output) {
  my_plot <- function(df, x, y, title, color) {
    ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
      geom_line(color = color, linewidth = 1) +
      labs(title = title) +
      theme_minimal()
  }
  sub_df <- reactive({
    y <- seq(0.001, 0.999, length.out = 300)
    dens <- dSUBXII(y, input$mu, input$sigma)
    cdf <- pSUBXII(y, input$mu, input$sigma)
    data.frame(y = y, dens = dens, cdf = cdf, surv = 1 - cdf, haz = dens/(1-cdf))
  })
  output$p1 <- renderPlot({ my_plot(sub_df(), "y", "dens", "PDF", "#3498db") })
  output$p2 <- renderPlot({ my_plot(sub_df(), "y", "cdf", "CDF", "#e74c3c") })
  output$p3 <- renderPlot({ my_plot(sub_df(), "y", "surv", "Survival", "#27ae60") })
  output$p4 <- renderPlot({ my_plot(sub_df(), "y", "haz", "Hazard", "#2c3e50") })
}
#-------------------------------------------------
shinyApp(ui, server)
#-------------------------------------------------