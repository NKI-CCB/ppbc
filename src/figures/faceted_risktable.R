#' Faceted risktables with ggplots
#'
#' @param ggkm A ggplot object generated with ggsurvplot with facet.by
#' @param br A vector of breaks for the Time axis
#' @param faceted.by Character, the response variable used with facet.by
#' @param text.size Numeric, the size of the risk table text
#' @param returndata Logical, whether to return the data plotted by the risk table
#'
#' @return A ggplot risktable with a layout comparable to ggkm
#' @export
#'
#' @examples
faceted_risktable <- function(ggkm, br=seq(0,300,50), faceted.by,
                              text.size = 12, returndata = F,
                              colors = c("NP-BC"="#521262", "Pr-BC"="#AA96DA",
                                         "PP-BCdl"="#112D4E", "PP-BCpw"="#46CDCF")
                              ){
  # Extract the table from the ggplot object
  df <- ggkm$data
  
  # Treating zero as 1 allows us to generate breaks starting from 0 without NAs
  df <- df %>%
    mutate(time1=ifelse(time == 0,1,time))
  
  # Create a series of time breaks
  df$breaks <- cut(df$time1, breaks=br)
  
  # Format the time breaks based on the start and convert to numeric
  df <- df %>%
    relocate(breaks, .after=time) %>%
    tidyr::separate(breaks, into = c("Time", NA), sep = ",", remove = F) %>%
    mutate(Time = as.numeric(gsub("\\(", "",Time)))
  
  # Summarize the max at risk patients for each time interval by strata
  df <- df %>%
    group_by(strata, Time) %>%
    mutate(N.risk = max(n.risk), .after=n.risk) %>%
    select(time:N.risk, strata:time1) %>%
    arrange(strata, Time)
  
  # Remove extraneous columns to plot a single unique value per interval and strata
  df <- df[,c("Time","N.risk", faceted.by, "study_group", "strata")] %>%
    distinct()
  
  # Hacky way of filling in the zeros for strata with less follow up
  df <- df %>%
    pivot_wider(names_from = Time, values_from = N.risk, values_fill = 0) %>%
    pivot_longer(cols = c(-study_group, -strata, -{{faceted.by}}),
                 names_to = "Time", values_to = "N.risk") %>%
    mutate(Time = as.numeric(Time))
  
  # Generate facet headings based on the facet group
  df$response <- paste(faceted.by, df[,faceted.by,drop=T])
  
  if(returndata){return(df)}
  
  # Generate the plot itself
  df %>%  
    ggplot(aes(x=Time, y=study_group, color = study_group)) + 
    geom_text(aes(label = N.risk), show.legend = F) +
    facet_wrap(~response) + theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill="white", size = 1),
          axis.text.x = element_text(color="black", size = text.size),
          axis.text.y = element_text(color="black", size = text.size),
          axis.title = element_text(color="black", size = text.size),
          strip.text = element_text(color="black", size = text.size))+
    ggtitle(paste("Risk table for", faceted.by)) +
    xlab("Time") + ylab("") +
    scale_color_manual(values = colors)
}