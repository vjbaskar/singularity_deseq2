getPrjName <- function(projectCSV){
    temp = strsplit(projectCSV, "[.]")[[1]]
    temp = temp[-length(temp)]
    temp
    pname = paste0(temp, sep="", collapse = ".")
    return(pname)
}

info <- function(text){
    text =  sapply(text, function(x) paste0(x, sep=" ", collapse = " "))
    d = format(Sys.time(), "%d/%m/%y %T")
    cat("[ Info: ", d, "] ")
    cat(text,"\n")
}

errorMessage <- function(text){
    text =  sapply(text, function(x) paste0(x, sep=" ", collapse = " "))
    d = date()
    cat(" [ Error: ", d, "]\n")
    cat(text,"\n")
    cat("Quitting")
    message(" [ Error: ", d, "]\n")
    message(text,"\n")
    cat("Quitting")
    quit(save="no", status=1)
}


createFileLink <- function(fileName){
    return(paste0(projectDir,"/", fileName))
}

stripNALines <- function(ss){
    # Remove NA cols
    temp = apply(ss, 2, unique)
    toRemoveCols = names(temp [ is.na(temp) ])
    if(length(toRemoveCols) > 0 ){
       # ss = dplyr::select(ss, ! ( colnames(ss) %in% toRemoveCols) )
       ss = ss[, ! ( colnames(ss) %in% toRemoveCols) ]
    }
    # Remove data if no runid or sampleid present
    ss = ss [ !is.na(ss[,2]), ]
    ss = ss [ !is.na(ss[,1]), ]
    return(ss)
}

cleanSheet <- function(ss){
    ss = stripNALines(ss)
    for(j in 1:ncol(ss)){
        ss[,j] <- as.character(ss[,j])
    }
    
    head(ss)
    for(i in 1:nrow(ss)){
        for(j in 1:ncol(ss)){
            t <- as.character(ss[i,j])
            keepOnlyAlphaNumeric(t)
            
            ss[i,j] <- keepOnlyAlphaNumeric(t)
        }
    }
    cn <- colnames(ss)
    cn <- sapply(cn, keepOnlyAlphaNumeric)
    cn <- toupper(cn)
    colnames(ss) <- cn
    return(ss)
}



keepOnlyAlphaNumeric <- function(text){
    specialChars <- "a1~!@#$%^&*(){}_+:\"<>?./;'[]-= "
    # text = gsub("[[ ].$]", "",text)
    text = gsub("[ ]", "_",text)
    text = gsub("[^[:alnum:]^,^_^-]", "",text)
    text = gsub("_*$", "",text)
    return(text)
}

# Check for reqd args.
checkvars <- function(reqdArgs){
    n=0
    for(i in reqdArgs){
        if(nchar(i)==0){
            n = n + 1 
        }
    }
    if(n > 0){
        errorMessage("Provide all the above args")
    }
}

Help <- function(helpVals){

    formatHelp <- function(temp, n){
    temp [ temp[,2] %in% n, 4 ] <- paste0( temp [ temp[,2] %in% n, 4 ], "\t", "<<--")
    
        for(i in 1:nrow(temp)){
        	
            cat(paste0("-",temp[i,1]),paste0("|--", temp[i,2]), temp[i,4],"\n")
        }
    }
    n = c()
    reqFields = helpVals[helpVals[,3] == 1,]
    reqdArgs = unique(reqFields[,2])
    for(i in reqdArgs){
        temp = tryCatch(
            { 
                Sys.getenv(i)
                temp = NA
            }, 
            error = function(x) {
                info(paste0("Var not Found: ", x))
                return(i)
            }
        )
        n = c(n, temp)
        #n = n + temp
    }
    n = unique(n)
    temp = Sys.getenv(reqdArgs, unset = NA)
    n = names(temp [ is.na(temp)])
   #	message(n)
    if(length(n) >= 1) {
        notgiven = as.data.frame(helpVals [ ! helpVals[,2] %in% n, ])
        given = helpVals [ helpVals[,2] %in% n, ]
        
        cat("Provide all the following options: \n")
        cat("---------------------------------- \n")
        formatHelp(reqFields, n)
 #       quit(save="no", status=2)
    }
}

