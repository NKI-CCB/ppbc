#' Faceted risktables with ggplots
#'
#' @param ggkm A ggplot object generated with ggsurvplot with facet.by
#' @param br A vector of breaks for the Time axis
#' @param faceted.by Character, the response variable used with facet.by
#' @param text.size Numeric, the size of the risk table text
#' @param returndata Logical, whether to return the data plotted by the risk table
#' @param returnMissing Logical, return a version of the data frame with missing intervals as NA (debug)
#' @param colors Character, a named vector of colors for the values within the table with names matching curve.var
#' (can be overwritten with scale_color_manual(values =c(rep("black", length(levels(curve.var))))) )
#'
#' @return A ggplot risktable with a layout comparable to ggkm
#' @export
#'
#' @examples faceted_risktable(IG_kaplan_PPBC_OS, faceted.by = "IG",curve.var = "study_group") +
#' scale_color_manual(values =c(rep("black", 3))) + theme_bw() +
#' theme(axis.text.y = element_text(colour = c(npbc = "#521262", prbc="#AA96DA",ppbcpw="#46CDCF")))
faceted_risktable <- function(ggkm, br=seq(0,300,50), faceted.by, curve.var="study_group",
                              text.size = 12, returndata = F, returnMissing=F,
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
  
  # Set aside a key-value pairing of strata and plotting variables
  plotvars <- df %>% ungroup() %>%
    select(strata, {{faceted.by}}, {{curve.var}}) %>%
    distinct()
  
  # Remove extraneous columns to plot a single unique value per interval and strata
  df <- df %>%
    select(Time, N.risk, strata) %>%
    distinct()
  
  # Fill in the time intervals with no events are omitted by cut() with NA
  df <- df %>%
    mutate(Time = as.factor(Time)) %>%
    tidyr::complete(Time, strata) %>%
    left_join(., plotvars, by = "strata") %>%
    arrange(strata, Time) %>% distinct() %>%
    mutate(Time = as.numeric(as.character(Time)))
  if(returnMissing){return(df)}
  
  # Replace the NAs with 0 if empty due to lack of follow up
  timeMax <- max(df$Time)
  df <- df  %>%
    mutate(N.risk = if_else(is.na(N.risk) & Time == timeMax, 0, N.risk)) %>%
    arrange(strata, Time)
  
  # Replace NA with the next interval if an earlier interval was skipped due to lack of events
  df$N.risk[is.na(df$N.risk)] <- df$N.risk[which(is.na(df$N.risk))+1]
  
  # Generate facet headings based on the facet group
  df$response <- df[,faceted.by,drop=T]
  
  if(returndata){return(df)}
  
  # Generate the plot itself
  df %>%  
    ggplot(aes(x=Time, y=get(curve.var), color = get(curve.var))) + 
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
    xlab("Time") + ylab(str_replace_all(curve.var, "_", " ")) +
    labs(str_replace_all(curve.var, "_", "")) +
    scale_color_manual(values = colors)
}