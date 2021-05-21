library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Training course example"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      selectInput("select", h3("Number of locations"), 
                  choices = list("1" = 1, "2" = 2,
                                 "5" = 3, "10"=4), selected = 1),
      
      selectInput("select2", h3("Number of months"), 
                  choices = list("1" = 1, "3" = 2,
                                 "6" = 3, "12"=4), selected = 1),
      
      selectInput("select_diagnosis", h3("Method of diagnosis"), 
                  choices = list("Suspected" = 1, "Confirmed (RDT)" = 2,
                                 "Confirmed (microscopy)" = 3), selected = 1),
      
      selectInput("select_haplotype", h3("Haplotype"), 
                  choices = list("loc1" = 1, "loc2" = 2), selected = 1),
      
      width=2
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot"),
      width=10
      
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  set.seed(2)
  
  N_locs <- 10
  locs <- cbind(runif(N_locs, 22, 24),
                runif(N_locs, -17, -15))
  
  #pick index to be centre
  centre_i <- 1
  #get distances of other points from this one
  loc_dists <- apply(locs, 1, function(row) sqrt(sum((row - locs[centre_i, ])^2)))
  #reorder locs
  locs <- locs[order(loc_dists), ]
  
  
  locs_plot <- list(1,
                    2:3,
                    c(1, 3, 5, 7, 9),
                    1:10)
  
  months_plot <- list(1, 
                      c(1, 5, 9),
                      c(1, 3, 5, 7, 9, 11),
                      1:12)
  
  plot(locs, asp=1)
  
  location_lists <- list()
  
  n_mutants_start <- c(100, 50, 50, 10, 10, 1, 1, rep(1, 3))
  relative_fitness <- c(1.05, 1.05, 1.02, 1.03, 1.01, rep(1, 5))
  
  n_generations <- 12
  
  
  
  for(i_loc in 1:N_locs){
    n_parasites_wildtype <- c()
    n_parasites_mutant <- c()
    n_parasites_wildtype[1] <- 1000
    n_parasites_mutant[1] <- n_mutants_start[i_loc]
    
    
    #distribution to replicate each generation
    #construct so expected value is 1
    f_finish_distribution <- function(offspring_distribution, expected_value=1){
      offspring_distribution[2] <- expected_value - sum(2:5 * offspring_distribution[3:6])
      offspring_distribution[1] <- 1 - sum(offspring_distribution[2:6])
      return(offspring_distribution)
    }
    
    offspring_distribution_wildtype <- rep(NA, 5)
    offspring_distribution_wildtype[3:6] <- c(0.1, 0.05, 0.05, 0.01)
    offspring_distribution_wildtype <- f_finish_distribution(offspring_distribution_wildtype)
    # print(offspring_distribution_wildtype)
    
    offspring_distribution_mutant <- rep(NA, 5)
    offspring_distribution_mutant[3:6] <- c(0.1, 0.05, 0.05, 0.01)
    offspring_distribution_mutant <- f_finish_distribution(offspring_distribution_mutant, relative_fitness[i_loc])
    # print(offspring_distribution_mutant)
    
    n_parasite_list <- list(wildtype=n_parasites_wildtype,
                            mutant=n_parasites_mutant)
    
    offspring_distribution_list <- list(wildtype=offspring_distribution_wildtype,
                                        mutant=offspring_distribution_mutant)
    
    for(i in 2:n_generations){
      for(j in 1:2){
        n_parasite_list[[j]][i] <- sum(sample(0:5, n_parasite_list[[j]][i-1], replace = TRUE, prob=offspring_distribution_list[[j]]))
        
      }
    }
    
    #remove every third obs
    for(i in 1:2){
      n_parasite_list[[i]] <- n_parasite_list[[i]][seq(1, n_generations, 1)]
    }
    
    # plot(n_parasite_list[[2]] / (n_parasite_list[[1]] + n_parasite_list[[2]]))
    location_lists[[i_loc]] <- n_parasite_list
  }
  
  n_months <- length(n_parasite_list[[1]])
  
  plot_df <- data.frame(x=rep(locs[, 1], each=n_months),
                        y=rep(locs[, 2], each=n_months),
                        loc=rep(1:N_locs, each=n_months),
                        n_wildtype=Reduce(c, lapply(location_lists, "[[", 1)),
                        n_mutant=Reduce(c, lapply(location_lists, "[[", 2)),
                        t=rep(1:n_months, N_locs))
  
  
  plot_df$prop_mutant <- plot_df$n_mutant / (plot_df$n_mutant + plot_df$n_wildtype)
  plot_df$prop_mutant <- plot_df$prop_mutant + rnorm(n=length(plot_df$prop_mutant), mean=0, sd=0.001)
  
  plot_df$loc <- as.factor(plot_df$loc)
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    plot_df_use <- plot_df
    which_use <- plot_df$loc %in% locs_plot[[as.numeric(input$select)]]
    plot_df_use$x[!which_use] <- NA
    plot_df_use$prop_mutant[!which_use] <- NA
    
    # map <- ggplot(plot_df[plot_df$loc %in% locs_plot[[as.numeric(input$select)]], ]) +
    # map <- ggplot(plot_df_use[plot_df_use$t %in% months_plot[[as.numeric(input$select2)]], ]) + 
    map <- ggplot() + 
      geom_point(data = plot_df, aes(x=x, y=y), col="grey", size=3) +
      geom_point(data = plot_df_use, aes(x=x, y=y, col=loc), size=3) +
      scale_color_discrete(drop=FALSE) + 
      theme_classic() + 
      theme(legend.position = "none") + 
      xlab("Longitude") + 
      ylab("Latitude") + 
      ylim(-17, -15) + 
      xlim(22, 24)
    
    # time_series <- ggplot(plot_df[plot_df$loc %in% locs_plot[[as.numeric(input$select)]], ]) + 
    time_series <- ggplot(plot_df_use[plot_df_use$t %in% months_plot[[as.numeric(input$select2)]], ]) + 
      geom_point(aes(x=t, y=prop_mutant, col=loc), size=2) +
      geom_line(aes(x=t, y=prop_mutant, col=loc)) + 
      scale_color_discrete(drop=FALSE) + 
      theme_linedraw() + 
      theme(legend.position = "none") +
      xlab("Month") + 
      ylab("Proportion of samples with mutation") + 
      ggtitle(NULL) + 
      ylim(-0.01, 0.2) + 
      scale_x_continuous(breaks=1:12) + 
      NULL
    
    out_plot <- patchwork::wrap_plots(time_series, map, widths=c(1, 0.75))
    
    return(out_plot)
    
  })
  
}


shinyApp(ui = ui, server = server)
