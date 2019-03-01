identifiers_for_topic_analysis <- list( 
  #BSH
  c(#"BSHINIT",
    #"Interpreter", 
     "isWrapperType", 
    "setNameSpace",
    "BeanShell",
    #"getResource", "importObject",
    #"pathToFile",  #"strictJava",
   "setNameSpace",  "eval",   "pwd" ), 
  
  #Regular Expression
  c("getMainRuleSet", #"expression", 
    "startRegexp", "endRegexp",  #"MATCH_TYPE_CONTEXT", 
   # "MATCH_TYPE_RULE",
    "setSearchMatcher" ,
    "pattern",
  "getLeadingWhiteSpace"), 
  
  
  #Text Area
  c( "lineHighlight" ,"highlight", "indentLine" ,"autoIndent" , 
    "getSelectedText","moveCaretPosition"),#"isJavaBaseAssignable", "resolveJavaMethod","resolveJavaField"  ) , #, 
  
  
  #XML support
  #"parseXML",  
  #"Registers" , "registerTransferableService", "saveRegisters", "setRegister",
  
  #User Interface
  c(
    "componentShown" ,"ScrollLayout",
    "focusedComponent",  "actionBar","JColorChooser" ,
    "showToolbars"),
  #"editPane", 
  #,"OPEN_DIALOG", "SAVE_DIALOG","addDockableWindow",
  
  #Core
  c("settingsDirectory",   #"zipFile", 
    "processKeyEvent", "getKeyEventInterceptor", 
    "createVFSSession", 
  "bufferLoaded",
     "addBufferListener"),#"MESSAGE",
  
  #plugins
  c("activatePlugin"  , "addPluginJAR",#"getJARCacheDirectory",
    "PluginOptions", 
    "pluginMgr","pluginDirectory" ,
   "getPluginVersion"),
  
  #Macros
  # "loadMacros",  "runScript", 
  
  #Search & Replace
  c("doBackwardSearch", "doForwardSearch", "searchStart",
    "searchField",  "hyperSearch" , "replaceSelection")# , 
)