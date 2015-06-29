# STATS GBMPRED extension command

#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2014
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# author=  'jkp, IBM'
# version=  '1.0.0'
# history
# 03-30-2013 original version

helptext = 'STATS GBMPRED MODELFILE=filespec
ID=varname
/SAVE DATASET=dataset name /INCLUDEIND = YES or NO
/OPTIONS BESTTREES=YES or NO NTREES=integer or list of integers
PREDSCALE = LINK or RESPONSE
/HELP.

This procedure calculates predicted values for new data using
results saved from the STATS GBM extension command.

All of the predictors from the original model must appear in the
dataset on which predictions will be made.  Missing values in
the predictors will be treated the same way as they were in
the estimation process.

MODELFILE specifies to get the mode from a file saved by
STATS GBM.  If it is not specified, the estimation results are
taken from the currently loaded R workspace.

The command output is a new dataset.  ID optionally specifies an
ID variable to facilitate merging this with the input data.

DATASET specifies a name for the predicted values dataset.  The dataset
name must not already be in use.  The dataset variables will be
the ID variable, if specified,
the predicted values
the independent variables if INCLUDEIND = YES.

The number of trees used for the predictions is controlled by
BESTTREES and NTREES.  If a best trees value was computed, it can
be used by specifying BESTTREES=YES.  Otherwise, NTREES specifies
the number of trees or a list of numbers of trees to use.  If both
BESTTREES and NTREES are specified, the besttrees value is prepended
to the NTREES list.  If neither is specified, the number of trees
specified at estimation time is used.
There will be one variable in the output dataset for each number of trees
specified.  If the dependent variable in the model is named
Y, the prediction variables will have names like Y_110, Y_500,...
where the suffix is the number of trees used for that set of predictions.

PREDSCALE specifies what is being predicted.  If the model has the
form Y= f(x), LINK creates predictions for f(x), i. e., the link
function.  For Bernoulli it calculates the log odds; for Poisson,
the log of counts, and for coxph, the log hazard.  
If RESPONSE is specified, predictions are converted
to the outcome scale.  This only affects  Bernoulli and Poisson.

/HELP displays this text and does nothing else.
'

library(gbm)

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}


doGbmPred <- function(dataset, modelfile=NULL, id=NULL,
    includeind=TRUE, besttrees=NULL, ntrees=NULL, predscale="link") {
    # create a new dataset containing predicted values and other relevant
    # values for variables in dataset
    # modelfile optionally specifies a file containing the workspace
    # id specifies an optional id variable
    # includeind includes the independent variables in the new dataset or not
    # besttree TRUE means use the previously calculated best number of trees
    # ntrees specifies a list of tree values to use
    # predscale specifies link or response to control what is predicted

    setuplocalization("STATS_GBMPRED")		

    tryCatch(library(gbm), error=function(e){
        stop(gtxtf("The R %s package is required but could not be loaded.","gbm"),call.=FALSE)
        }
    )
    if (is.null(modelfile))
        model = "workspace"
    if (is.null(ntrees)) {
        besttrees = TRUE}
    else {
        besttrees = FALSE
    }
    if (!is.null(modelfile))
        tryCatch(load(modelfile), error=function(e) {
            stop(sprintf(gtxt("Model file was not found or could not be read: %s"),modelfile), call.=FALSE)})
    # ensure that we have output from STATS GBM.  Only checking the result object
    tryCatch(class(res) == "gbm", error = function(e) {
        stop(gtxt("The specified file or workspace does not contain a gbm model estimated by STATS_GBM"),
        call.=FALSE)
        })
    if (besttrees) {
        besttree = tryCatch(modelproperties$bestiter, error = function(e) {
        stop(gtxt("Best number of trees value was not saved.  Use boostplot to save it."), 
        call.=FALSE)})
        ntrees = c(besttree, ntrees)   # if ntrees is NULL, it disappears from the vector
    }
    # If nothing was specified, use the estimation-time number of trees
    if (is.null(ntrees))
        ntrees = res$n.trees

    # display the previously estimated model
    StartProcedure(gtxt("Generalized Boosted Regression Predictions"), "STATSGBMPRED")
    spsspivottable.Display(settingsdf, 
        title = gtxt("Settings"),
        templateName = "GBMSUMMARY",
        outline=gtxt("Summary"),
        caption = gtxt("Results calculated by the R gbm procedure")
    )
    # get data for predictions from SPSS
    indep = res$var.names
    if (modelproperties["offset"] != "<None>")
        print(sprintf(gtxt("Note: predicted values do not include the offset variable: %s"),
            modelproperties["offset"]))
    # fetch data the same way as for estimation
    allvars = indep
    if (!is.null(id))
        allvars[length(allvars)+1] = id

    dta <- spssdata.GetDataFromSPSS(allvars, missingValueToNA = TRUE, 
        keepUserMissing=modelproperties["missingvalues"], factorMode = "levels")    
        
    predvalues = data.frame(predict(res, newdata=dta, n.trees=ntrees, type=predscale))
    
    spssdict = spssdictionary.GetDictionaryFromSPSS()
    # build output dataset dictionary
    # if including independent variables, reserve names in dictionary first
    slot = 1
    if (includeind) {
        nameset = indep
    } else {
        nameset = list()
    }
    dictlist = list()
    if (!is.null(id)) {
        dictlist[slot] = list(spssdict[match(id, spssdict["varName",])])
        nameset[length(nameset) + 1] = id
        slot = slot + 1
    }
    # Dependent variable not assumed to be in active dataset
    depvarspec = modelproperties["depvar"] # list of properties
    # Fill in entries for each predicted variable
    # Names are generated as estimation dep variable with a tree count suffix
    depvarinputspec = depvarspec
    for (i in 1:length(ntrees)) {
        depvarname = getname(depvarinputspec, ntrees[i], nameset)
        depvarspec[[1]][1] = depvarname   # replace name in model with suffixed name
        dictlist[slot] = depvarspec
        nameset[slot] = depvarname
        slot = slot + 1
    }

    if (includeind) { # assuming we can use current independent variable properties
        for (v in indep) {
            dictlist[slot] = list(spssdict[match(v, spssdict["varName",])])
            slot = slot + 1
        }
    }
    spsspkg.EndProcedure()
    dictlist = do.call(spssdictionary.CreateSPSSDictionary, dictlist)
    tryCatch(spssdictionary.SetDictionaryToSPSS(dataset, dictlist),
        error = function(e) {spssdictionary.EndDataStep(); stop(as.character(e), call.=FALSE)})
    # build results data frame
    if (!is.null(id))
        predvalues = data.frame(dta[id],predvalues)
    if (includeind) {
        if (is.null(id)) {
            predvalues = data.frame(predvalues, dta)
        } else {
        predvalues = data.frame(predvalues, dta[-length(dta)])
        }
    }
    spssdata.SetDataToSPSS(dataset, predvalues)
    spssdictionary.EndDataStep()
}

getname = function(depvarspec, ntrees, nameset) {
    root = depvarspec[[1]][1]
    roottrial = paste(root, ntrees, sep="_")
    trial = 1
    while (roottrial %in% nameset) {
        roottrial = paste(root, ntrees, trial, sep="_")
        trial = trial + 1
    }
    return(roottrial)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

gtxt <- function(...) {
    return(gettext(...,domain="STATS_GBMPRED"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_GBMPRED"))
}

Run <- function(args) {
#Execute the STATS GBMPRED extension command
    
    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
    spsspkg.Template("MODELFILE", subc="", ktype="literal", var="modelfile"),
    spsspkg.Template("ID", subc="", ktype="existingvarlist", var="id"),
    spsspkg.Template("BESTTREES", subc="OPTIONS", ktype="bool", var="besttrees"),
    spsspkg.Template("NTREES", subc="OPTIONS", ktype="int", var="ntrees", islist=TRUE),
    spsspkg.Template("PREDSCALE", subc="OPTIONS", ktype="str", var="predscale", 
        vallist=list("link", "response")),
    spsspkg.Template("DATASET", subc="SAVE", ktype="varname", var="dataset"),
    spsspkg.Template("INCLUDEIND", subc="SAVE", ktype="bool", var="includeind"),
    spsspkg.Template("HELP", subc="", ktype="bool")
    ))


# A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doGbmPred")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}