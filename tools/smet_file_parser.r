#####################################################################################
#  Copyright 2019 Michael Reisecker                                                 #
#####################################################################################
#  This is free software: you can redistribute and/or modify it under the terms of  #
#  the GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
#####################################################################################

#SYNOPSIS: smet <- read_smet("<smetfile>"); print(smet@data$TA[[1]])

#Methods to parse a SMET file. Note that there exists a more powerful and complicated
#R package for this called RSMET (CRAN archive).

library(methods) #for setClass

########################################
#  SMET CLASS DEFINITION               #
########################################

setClass("SMET", slots = list(header = "list", data = "data.frame", time = "POSIXlt", fields = "character",
    description = "data.frame", color = "data.frame", station_id = "character"))

invisible( #constructor: read input file
setMethod("initialize", signature = "SMET",
	definition = function(.Object, infile) {
		.Object <- read.smet(.Object, infile)
	  	return(.Object)
	}
))

invisible( #no output on definition
setGeneric("read.smet",
	function(smet, infile) {
		standardGeneric("read.smet")
	}
))

invisible(
setMethod("read.smet", signature = "SMET",
	function(smet, infile) { #read and process a smet file
		raw = get_raw_smet(infile)

		smet@header = raw[[1]]
		smet@data = as.data.frame(raw[[2]])
		smet@station_id <- get_header_value(smet, "station_id")
		smet@fields <- get_header_value(smet, "fields")

		colnames(smet@data) <- smet@fields

		#convert time strings to time list and save to dedicated slot:
		smet@time <- as.POSIXlt(smet@data$time, format = "%Y-%m-%dT%H:%M:%OS")
		smet@data$time <- NULL #remove time from dataframe
		return(smet)
	}
))

########################################
#  SMET METHODS                        #
########################################

invisible(
setGeneric("convert.units",
	function(smet) {
		standardGeneric("convert.units")
	}
))

invisible(
setMethod("convert.units", signature = "SMET", #bind this function to the SMET class
	function(smet) { #convert SI to human friendly
		if ("TA" %in% names(smet@data)) smet@data$TA = k2c(smet@data$TA) #air temperature
		if ("TA_MAX" %in% names(smet@data)) smet@data$TA_MAX = k2c(smet@data$TA_MAX)
		if ("TA_MIN" %in% names(smet@data)) smet@data$TA_MIN = k2c(smet@data$TA_MIN)
		if ("TD" %in% names(smet@data)) smet@data$TD = k2c(smet@data$TD) #dew point
		if ("HI" %in% names(smet@data)) smet@data$HI = k2c(smet@data$HI) #heat index
		if ("WB" %in% names(smet@data)) smet@data$WB = k2c(smet@data$WB) #wet bulb
		if ("WC" %in% names(smet@data)) smet@data$WC = k2c(smet@data$WC) #wind chill
		if ("RH" %in% names(smet@data)) smet@data$RH = ratio2percent(smet@data$RH) #humidity
		if ("RH_MAX" %in% names(smet@data)) smet@data$RH_MAX = ratio2percent(smet@data$RH_MAX)
		if ("RH_MIN" %in% names(smet@data)) smet@data$RH_MIN = ratio2percent(smet@data$RH_MIN)
		if ("VW" %in% names(smet@data)) smet@data$VW = mps2kmph(smet@data$VW) #wind speed
		if ("VW_MAX" %in% names(smet@data)) smet@data$VW_MAX = mps2kmph(smet@data$VW_MAX) #wind gusts
		return(smet)
  }
))

invisible(
setGeneric("get.timeframe", #timeframe of dataset in days
	function(smet) {
		standardGeneric("get.timeframe")
	}
))

invisible(
setMethod("get.timeframe", signature = "SMET",
	function(smet) { #return difference in days:
		return(difftime(smet@time[length(smet@time)], smet@time[1], units = "days"))
  }
))

########################################
#  HELPER FUNCTIONS                    #
########################################

check_param <- function(smet, param) { #can the parameter be worked with at all?
	if (!param %in% names(smet@data)) {
		return(FALSE)
	}
	if (all( is.na(smet@data[[param]]) )) {
		return(FALSE)
	}
	if (all( smet@data[[param]] == 0, na.rm = TRUE )) { #could be good data (e.g. no precipitation)
		return(FALSE) #but the plots are still of little use
	}
	return(TRUE)
}

k2c <- function(xx) {
	return(xx - 273.15)
}

ratio2percent <- function(xx) {
	return(xx * 100)
}

mps2kmph <- function(xx) {
	return(xx * 3.6)
}

########################################
#  FILE READING                        #
########################################

get_raw_smet <- function(infile) { #read header and data lines from file system
	idx.datasection <- grep("\\[DATA\\]", readLines(infile), value = FALSE) #return index of line "[DATA]"
	data <- read.csv(file = infile, sep = "", skip = idx.datasection, head = FALSE, na.strings = "-999",
	    stringsAsFactors = FALSE, strip.white = TRUE)

	idx.headersection <- grep("\\[HEADER\\]", readLines(infile), value = FALSE) #return index of line "[HEADER]"
	header <- read.table(file = infile, sep = "=", skip = idx.headersection, nrows = idx.datasection - idx.headersection - 1,
	    stringsAsFactors = FALSE, strip.white = TRUE)

	return(list(header, data))
}

get_header_value <- function(smet, key) { #find a key in the header and return its value ("key = value")
	idx <- grep(key, smet@header[[1]])
	str <- smet@header[[2]][idx]
	return(strsplit(str, " ", fixed = TRUE)[[1]]) #unlist
}

########################################
#  WRAPPER                             #
########################################

read_smet <- function(inputfile, convert_units = TRUE) { #wrapper for class initialization
	smet <- new("SMET", infile = inputfile)
	if (convert_units) {
		smet <- convert.units(smet)
	}
	return(smet)
}

