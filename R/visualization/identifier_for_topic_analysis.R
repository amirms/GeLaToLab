identifiers_for_topic_analysis <- list( 
  #BSH
  c(#"BSHINIT",
    "Interpreter", #"invoke", #"isWrapperType", 
    #"setNameSpace","setAccessibility",
    # "setStrictJava",  "BSHForStatement","getMethodNames","getSourceFileInfo",
    #"BeanShell",
    #"getResource", "importObject",
    "pathToFile",  "strictJava", "setNameSpace", "exec", "eval",   "pwd",  "debug", 
     "BSHPackageDeclaration" ,"BSHBinaryExpression"),# "searchJarForClasses" ), 
  
  #Regular Expression
  c(#"escapeRule" , #"getMainRuleSet", "expression", "terminateChar", 
    "startRegexp", "endRegexp", "regexp",  "MATCH_TYPE_CONTEXT", "MATCH_TYPE_RULE",
    # "setSearchMatcher" ,
    "pattern", "matchType") , #"MATCH_TYPE_CONTEXT", "MATCH_TYPE_RULE"  , "terminateChar", 
  #"getLeadingWhiteSpace", #"match","escapeRegexp", 
  
  
  #Text Area
  c("caretLine", "findMatchingBracket" ,"highlight", "wrap", #"selectToMatchingBracket", #"lineHighlight" ,"highlight", "indentLine" ,"autoIndent" , 
    "getLineCount" ,  "caret", "getCaretPosition", #"getCaretLine",
    "getSelectedText"),#"isJavaBaseAssignable", "resolveJavaMethod","resolveJavaField"  ) , #"moveCaretPosition", 
  
  
  #XML support
  #"parseXML",  
  #"Registers" , "registerTransferableService", "saveRegisters", "setRegister",
  
  #User Interface
  #"preferredLayoutSize" , 
  c("toolbar", "menubar" ,
    "showPopupMenu", "componentShown" ,# "ScrollLayout",
    #"needFullRepaint", 
    "focusedComponent",  "actionBar",#"JColorChooser" ,
    "showToolbars"),
  #"editPane", 
  #,"OPEN_DIALOG", "SAVE_DIALOG","addDockableWindow",
  
  #Core
  c("settingsDirectory",   "zipFile", "processKeyEvent", "getKeyEventInterceptor", #"handleMessage", 
    # "fileVFS", #"createVFSSession", , #"sourceFile" , "queueAWTRunner",
    "invokeAction", "buffer",  "bufferCount"),# "bufferLoaded"),# "ZipInputStream"), "handleMessage",
  #"bufferCount", #"bufferLoaded", "bufferOpened",  #"queueAWTRunner", "createVFSSession",
  # "updateBufferStatus", 
  #   "addBufferListener", #"MESSAGE", "performOperationsInAWTThread",
  
  #plugins
  c("activatePlugin"  , #"addPluginJAR","getJARCacheDirectory", #"PluginOptions", 
    # "download", 
   # "pluginMgr",#"pluginDirectory" ,
    "author", "description", "version", "pluginSet", "getPluginJAR", "getPluginVersion","PluginOptions"),
  
  #Macros
  # "loadMacros",  "runScript", 
  
  #Search & Replace
  c("doBackwardSearch", "doForwardSearch", "searchStart",
    "searchField",  "replace" , "replaceSelection")#"hyperSearch" , 
)