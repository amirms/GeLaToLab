# Used and Modified from code written by 
# http://www.stat.cmu.edu/~cshalizi/350/hw/01/01.R

# Read in a source text file
# Input: filename
# Calls: XML package (available from CRAN)
#        strip.text()
# Output: vector of character strings, giving the words in order
read.doc <- function(filename, lengthLowerBound=4) {
  
  fulltext <- readLines(filename, encoding="UTF-8")
  # ASSUMES: this should be a SINGLE character string
  #print("reached here")
  print(filename)
  text <- strip.java.text(fulltext, lengthLowerBound) # Turn into a vector of strings
  return(text)
}

read.text <- function(filename) {
  
  fulltext <- readChar(filename, file.info(filename)$size)  

#   fulltext <- gsub("[\n\t]+"," ", fulltext)
  
  return(fulltext)
}

# Feel free to use this IF you comment it
read.directory <- function(dirname, pattern="*", verbose=FALSE) {
  srcfiles = list()
  filenames = list.files(path = dirname, pattern = pattern, all.files = FALSE,
                         full.names = TRUE, recursive = TRUE,
                         ignore.case = TRUE, include.dirs = TRUE, no.. = TRUE)
#   srcfiles.names = c(0)
  for (i in 1:length(filenames)) {
    if(verbose) {
      print(filenames[i])
    }
    #srcfiles.names <- cbind(, srcfiles.names)
    srcfiles[[filenames[i]]] = read.doc(filenames[i])
    
  }
  
  #names(srcfiles) <- srcfiles.names
  
  return(srcfiles)
}

read.text.directory <- function(dirname, pattern="*", verbose=FALSE) {
  srcfiles = list()
  filenames = list.files(path = dirname, pattern = pattern, all.files = FALSE,
                         full.names = TRUE, recursive = TRUE,
                         ignore.case = TRUE, include.dirs = TRUE, no.. = TRUE)
#   srcfiles.names = c()
  for (i in 1:length(filenames)) {
    if(verbose) {
      print(filenames[i])
    }
#     srcfiles.names <- cbind(, basename(filenames[i]))
    srcfiles[[filenames[i]]] = read.text(filenames[i])
    
  }
  
#   names(srcfiles) <- srcfiles.names
  
  return(srcfiles)
}



tryTolower = function(x)
{
  # create missing value
  # this is where the returned value will be
  y = NA
  # tryCatch error
  try_error = tryCatch(tolower(x), error = function(e) e)
  # if not an error
  if (!inherits(try_error, "error"))
    y = tolower(x)
  return(y)
}

# Turn a string into a vector of words
# for comparability across bags of words, also strip out punctuation and
# numbers, and shift all letters into lower case
# Input: character string
# Output: vector of words (charaacter strings)
strip.text <- function(txt) {  
  
  #txt <- tolower(txt)
  
  #change other non-alphanumeric characters to spaces
  txt <- gsub("[^a-zA-Z0-9]"," ",txt)
  
  # change digits to #
  txt <- gsub("[0-9]+","#",txt)
  
  # split and make one vector
  
  txt <- strsplit(txt," ")
  
  #
  pat <- "(?<=[[:lower:]])(?=[[:upper:]])"
  txt <- lapply(txt, strsplit, pat, perl=TRUE)
  
  txt <- unlist(txt)
  
  # convert to lowercase
  txt <- tryTolower(txt)
    
  # remove words of one character or less
  txt <- txt[nchar(txt) > 1]
  
  # remove numeric data
  txt <- txt[!is.numeric(txt)]
  return(txt)
}



pli.reserved.words <- c("ALL", "ALTER", "AND", "ANY", "ARRAY", "ARROW", "AS", "ASC", "AT",
                        "BEGIN", "BETWEEN", "BY", 
                        "CASE", "CHECK", "CLUSTERS", "CLUSTER", "COLAUTH", "COLUMNS", "COMPRESS", "CONNECT", "CRASH", "CREATE", "CURRENT",
                        "DECIMAL", "DECLARE", "DEFAULT", "DELETE", "DESC", "DISTINCT", "DROP",
                        "ELSE", "END", "EXCEPTION", "EXCLUSIVE", "EXISTS",
                        "FETCH", "FORM", "FOR", "FROM",
                        "GOTO", "GRANT", "GROUP",
                        "HAVING",
                        "IDENTIFIED", "IF", "IN", "INDEXES", "INDEX", "INSERT", "INTERSECT", "INTO", "IS",
                        "LIKE", "LOCK",
                        "MINUS", "MODE",
                        "NOCOMPRESS", "NOT", "NOWAIT", "NULL",
                        "OF", "ON", "OPTION", "OR", "ORDER", "OVERLAPS",
                        "PRIOR", "PROCEDURE", "PUBLIC",
                        "RANGE", "RECORD", "RESOURCE", "REVOKE",
                        "SELECT", "SHARE", "SIZE", "SQL", "START", "SUBTYPE",
                        "TABAUTH", "TABLE", "THEN", "TO", "TYPE",
                        "UNION", "UNIQUE", "UPDATE", "USE",
                        "VALUES", "VIEW", "VIEWS",
                        "WHEN", "WHERE", "WITH")

pli.keywords <- c("ADD", "AGENT", "AGGREGATE", "ARRAY", "ATTRIBUTE", "AUTHID", "AVG",
                  "BFILE_BASE", "BINARY", "BLOB_BASE", "BLOCK", "BODY", "BOTH", "BOUND", "BULK", "BYTE",
                  "C", "CALL", "CALLING", "CASCADE", "CHAR", "CHAR_BASE", "CHARACTER", "CHARSETFORM", "CHARSETID", "CHARSET",
                  "CLOB_BASE", "CLOSE", "COLLECT", "COMMENT", "COMMIT", "COMMITTED", "COMPILED", "CONSTANT", "CONSTRUCTOR", 
                  "CONTEXT", "CONVERT", "COUNT", "CURSOR", "CUSTOMDATUM",
                  "DANGLING", "DATA", "DATE", "DATE_BASE", "DAY", "DEFINE", "DETERMINISTIC", "DOUBLE", "DURATION",
                  "ELEMENT", "ELSIF", "EMPTY", "ESCAPE", "EXCEPT", "EXCEPTIONS", "EXECUTE", "EXIT", "EXTERNAL",
                  "FINAL", "FIXED", "FLOAT", "FORALL", "FORCE", "FUNCTION",
                  "GENERAL",
                  "HASH", "HEAP", "HIDDEN", "HOUR",
                  "IMMEDIATE", "INCLUDING", "INDICATOR", "INDICES", "INFINITE", "INSTANTIABLE", "INT", "INTERFACE", 
                  "INTERVAL", "INVALIDATE", "ISOLATION",
                  "JAVA",
                  "LANGUAGE", "LARGE", "LEADING", "LENGTH", "LEVEL", "LIBRARY", "LIKE2", "LIKE4",
                  "LIKEC", "LIMIT", "LIMITED", "LOCAL", "LONG", "LOOP",
                  "MAP", "MAX", "MAXLEN", "MEMBER", "MERGE", "MIN", "MINUTE", "MOD", "MODIFY", "MONTH", "MULTISET",
                  "NAME", "NAN", "NATIONAL", "NATIVE", "NCHAR", "NEW", "NOCOPY", "NUMBER_BASE",
                  "OBJECT", "OCICOLL", "OCIDATETIME", "OCIDATE", "OCIDURATION", "OCIINTERVAL", "OCILOBLOCATOR", 
                  "OCINUMBER", "OCIRAW", "OCIREFCURSOR", "OCIREF", "OCIROWID", "OCISTRING", "OCITYPE", "ONLY", "OPAQUE", 
                  "OPEN", "OPERATOR", "ORACLE", "ORADATA", "ORGANIZATION", "ORLANY", "ORLVARY", "OTHERS", "OUT", "OVERRIDING",
                  "PACKAGE", "PARALLEL_ENABLE", "PARAMETER", "PARAMETERS", "PARTITION", 
                  "PASCAL", "PIPE", "PIPELINED", "PRAGMA", "PRECISION", "PRIVATE",
                  "RAISE", "RANGE", "RAW", "READ", "RECORD", "REF", "REFERENCE", "REM", "REMAINDER", 
                  "RENAME", "RESULT", "RETURN", "RETURNING", "REVERSE", "ROLLBACK", "ROW",
                  "SAMPLE", "SAVE", "SAVEPOINT", "SB1", "SB2", "SB4", "SECOND", "SEGMENT", "SELF", "SEPARATE", 
                  "SEQUENCE", "SERIALIZABLE", "SET", "SHORT", "SIZE_T", "SOME", "SPARSE", "SQLCODE", "SQLDATA", 
                  "SQLNAME", "SQLSTATE", "STANDARD", "STATIC", "STDDEV", "STORED", "STRING", "STRUCT", "STYLE", 
                  "SUBMULTISET", "SUBPARTITION", "SUBSTITUTABLE", "SUBTYPE", "SUM", "SYNONYM",
                  "TDO", "THE", "TIME", "TIMESTAMP", "TIMEZONE_ABBR", "TIMEZONE_HOUR", "TIMEZONE_MINUTE", "TIMEZONE_REGION", 
                  "TRAILING", "TRANSAC", "TRANSACTIONAL", "TRUSTED", "TYPE",
                  "UB1", "UB2", "UB4", "UNDER", "UNSIGNED", "UNTRUSTED", "USE", "USING",
                  "VALIST", "VALUE", "VARIABLE", "VARIANCE", "VARRAY", "VARYING", "VOID",
                  "WHILE", "WORK", "WRAPPED", "WRITE",
                  "YEAR",
                  "ZONE")


cobol.reserved.words <- c("ACCEPT", "ACCESS", "ADD", "ADDRESS", "ADVANCING", "AFTER", "ALL",
                          "ALLOWING", "ALPHABET", "ALPHABETIC", "ALPHABETIC--LOWER", "ALPHABETIC--UPPER", "ALPHANUMERIC",
                          "ALPHANUMERIC--EDITED", "ALSO", "ALTER", "ALTERNATE", "AND", "ANY", "APPLY", "ARE",  "AREA",
                          "AREAS", "ASCENDING", "ASSIGN", "AT", "AUTHOR", "AUTO", "AUTOMATIC", "AUTOTERMINATE",
                          "BACKGROUND-COLOR", "BATCH", "BEFORE", "BEGINNING", "BELL", "BINARY", "BINARY-CHAR", 
                          "BINARY-DOUBLE", "BINARY-LONG", "BINARY-SHORT", "BIT", "BITS", "BLANK", "BLINK", "BLINKINK",
                          "BLOCK", "BOLD", "BOOLEAN", "BOTTOM", "BY", "CALL", "CANCEL", "CD", "CF", "CH", "CHANGED",
                          "CHARACTER", "CHARACTERS", "CLASS", "CLOCK-UNITS", "CLOSE", "COBOL", "CODE", "CODE-SET",
                          "COL", "COLLATING", "COLUMN", "COMMA", "COMMIT", "COMMON", "COMMUNICATION", "COMP", "COMP-1",
                          "COMP-2", "COMP-3", "COMP-4", "COMP-5", "COMP-6", "COMP-X", "COMPUTATIONAL", "COMPUTATIONAL-1",
                          "COMPUTATIONAL-2", "COMPUTATIONAL-3", "COMPUTATIONAL-4", "COMPUTATIONAL-5", "COMPUTATIONAL-6",
                          "COMPUTATIONAL-X", "COMPUTE", "CONCURRENT", "CONFIGURATION", "CONNECT", "CONTAIN", "CONTAINS",
                          "CONTENT", "CONTINUE", "CONTROL", "CONTROLS", "CONVERSION", "CONVERTING", "COPY", "CORE-INDEX",
                          "CORR", "CORRESPONDING", "COUNT", "CRT", "CURRENCY", "CURRENT", "CURSOR", "DATA", "DATE", 
                          "DATE-COMPILED", "DATE-WRITTEN", "DAY", "DAY-OF-WEEK", "DB", "DB-ACCESS-CONTROL-KEY", 
                          "DB-CONDITION", "DB-CURRENT-RECORD-ID", "DB-CURRENT-RECORD-NAME", "DB-EXCEPTION", "DB-KEY",
                          "DB-RECORD-NAME", "DB-SET-NAME", "DB-STATUS", "DB-UWA", "DBCS", "DBKEY", "DE", "DEBUG-CONTENTS",
                          "DEBUG-ITEM", "DEBUG-LENGTH", "DEBUG-LINE", "DEBUG-NAME", "DEBUG-NUMERIC-CONTENTS", "DEBUG-SIZE",
                          "DEBUG-START", "DEBUG-SUB", "DEBUG-SUB-1", "DEBUG-SUB-2", "DEBUG-SUB-3", "DEBUG-SUB-ITEM",
                          "DEBUG-SUB-N", "DEBUG-SUB-NUM", "DEBUGGING", "DECIMAL-POINT", "DECLARATIVES", "DEFAULT",
                          "DELETE", "DELIMITED", "DELIMITER", "DEPENDENCY", "DEPENDING", "DESCENDING", "DESCRIPTOR",
                          "DESTINATION", "DETAIL", "DICTIONARY", "DISABLE", "DISCONNECT", "DISP", "DISPLAY", "DISPLAY-1",
                          "DISPLAY-6", "DISPLAY-7", "DISPLAY-9", "DIVIDE", "DIVISION", "DOES", "DOWN", "DUPLICATE",
                          "DUPLICATES", "ECHO", "EDITING", "EGI", "EJECT", "ELSE", "EMI", "EMPTY", "ENABLE", "END",
                          "END-ACCEPT", "END-ADD", "END-CALL", "END-COMMIT", "END-COMPUTE", "END-CONNECT", "END-DELETE",
                          "END-DISCONNECT", "END-DIVIDE", "END-ERASE", "END-EVALUATE", "END-FETCH", "END-FIND", "END-FINISH",
                          "END-FREE", "END-GET", "END-IF", "END-KEEP", "END-MODIFY", "END-MULTIPLY", "END-OF-PAGE", 
                          "END-PERFORM", "END-READ", "END-READY", "END-RECEIVE", "END-RECONNECT", "END-RETURN", "END-REWRITE",
                          "END-ROLLBACK", "END-SEARCH", "END-START", "END-STORE", "END-STRING", "END-SUBTRACT", "END-UNSTRING",
                          "END-WRITE", "ENDING", "ENTER", "ENTRY", "ENVIRONMENT", "EOL", "EOP", "EOS", "EQUAL", "EQUALS",
                          "ERASE", "ERROR", "ESI", "EVALUATE", "EVERY", "EXAMINE", "EXCEEDS", "EXCEPTION", "EXCLUSIVE",
                          "EXHIBIT", "EXIT", "EXOR", "EXTEND", "EXTERNAL", "FAILURE", "FALSE", "FD", "FETCH", "FILE",
                          "FILE-CONTROL", "FILLER", "FINAL", "FIND", "FINISH", "FIRST", "FLOAT-EXTENDED", "FLOAT-LONG",
                          "FLOAT-SHORT", "FOOTING", "FOR", "FOREGROUND-COLOR", "FREE", "FROM", "FULL", "FUNCTION", 
                          "GENERATE", "GET", "GIVING", "GLOBAL", "GO", "GOBACK", "GREATER", "GROUP", "HEADING",
                          "HIGH-VALUE", "HIGH-VALUES", "HIGHLIGHT", "I-O", "I-O-CONTROL", "ID", "IDENT", "IDENTIFICATION",
                          "IF", "IN", "INCLUDING", "INDEX", "INDEXED", "INDICATE", "INITIAL", "INITIALIZE", "INITIATE",
                          "INPUT", "INPUT-OUTPUT", "INSPECT", "INSTALLATION", "INTO", "INVALID", "IS", "JUST", "JUSTIFIED",
                          "KANJI", "KEEP", "KEY", "LABEL", "LAST", "LD", "LEADING", "LEFT", "LENGTH", "LESS", "LIMIT",
                          "LIMITS", "LINAGE", "LINAGE-COUNTER", "LINE", "LINE-COUNTER", "LINES", "LINKAGE", "LOCALLY",
                          "LOCK", "LOCK-HOLDING", "LOW-VALUE", "LOW-VALUES", "LOWLIGHT", "MANUAL", "MATCH", "MATCHES",
                          "MEMBER", "MEMBERSHIP", "MEMORY", "MERGE", "MESSAGE", "MODE", "MODIFY", "MODULES", "MOVE",
                          "MULTIPLE", "MULTIPLY", "NAMED", "NATIVE", "NEGATIVE", "NEXT", "NO", "NON-NULL", "NOT",
                          "NOTE", "NULL", "NUMBER", "NUMERIC", "NUMERIC-EDITED", "OBJECT-COMPUTER", "OCCURS", "OF",
                          "OFF", "OFFSET", "OMITTED", "ON", "ONLY", "OPEN", "OPTIONAL", "OPTIONS", "OR", "ORDER",
                          "OTHERWISE", "PACKED-DECIMAL", "PADDING", "PAGE", "PAGE-COUNTER", "PASSWORD", "PERFORM", "PF",
                          "PH", "PIC", "PICTURE", "PLUS", "POINTER", "POSITION", "POSITIONING", "POSITIVE", "PREVIOUS",
                          "PRINTING", "PRIOR", "PROCEDURE", "PROCEDURES", "PROCEED", "PROGRAM", "PROGRAM-ID", "PROTECTED",
                          "PURGE", "QUEUE", "QUOTE", "QUOTES", "RANDOM", "RD", "READ", "READERS", "READY", "REALM",
                          "REALMS", "RECEIVE", "RECONNECT", "RECORD", "RECORD-NAME", "RECORD-OVERFLOW", "RECORDING",
                          "RECORDS", "REDEFINES", "REEL", "REFERENCE", "REFERENCE-MODIFIER", "REFERENCES", "REGARDLESS",
                          "RELATIVE", "RELEASE", "RELOAD", "REMAINDER", "REMARKS", "REMOVAL", "RENAMES", "REORG-CRITERIA",
                          "REPLACE", "REPLACING", "REPORT", "REPORTING", "REPORTS", "REQUIRED", "RERUN", "RESERVE", 
                          "RESET", "RETAINING", "RETRIEVAL", "RETURN", "RETURN-CODE", "RETURNING", "REVERSE-VIDEO",
                          "REVERSED", "REWIND", "REWRITE", "RF", "RH", "RIGHT", "RMS-CURRENT-FILENAME", "RMS-CURRENT-STS",
                          "RMS-CURRENT-STV", "RMS-FILENAME", "RMS-STS", "RMS-STV", "ROLLBACK", "ROUNDED", "RUN", "SAME",
                          "SCREEN", "SD", "SEARCH", "SECTION", "SECURE", "SECURITY", "SEGMENT", "SEGMENT-LIMIT", "SELECT",
                          "SEND", "SENTENCE", "SEPARATE", "SEQUENCE", "SEQUENCE-NUMBER", "SEQUENTIAL", "SERVICE", "SET",
                          "SETS", "SIGN", "SIGNED", "SIZE", "SKIP1", "SKIP2", "SKIP3", "SORT", "SORT-MERGE", "SOURCE",
                          "SOURCE-COMPUTER", "SPACE", "SPACES", "SPECIAL-NAMES", "STANDARD", "STANDARD-1", "STANDARD-2",
                          "START", "STATUS", "STOP", "STORE", "STREAM", "STRING", "SUB-QUEUE-1", "SUB-QUEUE-2", 
                          "SUB-QUEUE-3", "SUB-SCHEMA", "SUBTRACT", "SUCCESS", "SUM", "SUPPRESS", "SYMBOLIC", "SYNC",
                          "SYNCHRONIZED", "TABLE", "TALLYING", "TAPE", "TENANT", "TERMINAL", "TERMINATE", "TEST",
                          "TEXT", "THAN", "THEN", "THROUGH", "THRU", "TIME", "TIMES", "TO", "TOP", "TRACE", "TRAILING",
                          "TRANSFORM", "TRUE", "TYPE", "UNDERLINE", "UNDERLINED", "UNEQUAL", "UNIT", "UNLOCK", "UNSIGNED", 
                          "UNSTRING", "UNTIL", "UP", "UPDATE", "UPDATERS", "UPON", "USAGE", "USAGE-MODE", "USE", "USING",
                          "VALUE", "VALUES", "VARYING", "VFU-CHANNEL", "WAIT", "WHEN", "WHERE", "WITH", "WITHIN", "WORDS",
                          "WORKING-STORAGE", "WRITE", "WRITERS", "ZERO", "ZEROES", "ZEROS")


easytrieve.reserved.words <- c("AFTER-BREAK", "AFTER-LINE", "AFTER-SCREEN", "AIM", "AND", "ATTR", 
                               "BEFORE", "BEFORE-BREAK", "BEFORE-LINE", "BEFORE-SCREEN", "BUSHU", "BY",
                               "CALL", "CASE", "CHECKPOINT", "CHKP", "CHKP-STATUS", "CLEAR", "CLOSE", "COL", 
                               "COLOR", "COMMIT", "CONTROL", "COPY", "CURSOR",
                               "D", "DECLARE", "DEFAULT", "DEFINE", "DELETE", "DENWA", "DISPLAY", "DLI", "DO", "DUPLICATE",
                               "E", "ELSE", "ELSE-IF", "END", "END-CASE", "END-DO", "END-IF", "END-PROC", "ENDPAGE", 
                               "ENDTABLE", "ENTER", "EOF", "EQ", "ERROR", "EXIT", "EXTERNAL", "EZLIB",
                               paste("F",  1:36, sep = ""), "FETCH", "FILE", "FILE-STATUS", "FILL", "FINAL", "FIRST", "FIRST-DUP", "FOR",
                               "GE", "GET", "GO", "GOTO", "GQ", "GR", "GT",
                               "HEADING", "HEX", "HIGH-VALUES",
                               "IDD", "IDMS", "IF", "IN", "INSERT",
                               "JOB", "JUSTIFY",
                               "KANJI-DATE", "KANJI-DATE-LONG", "KANJI-TIME", "KEY", "KEY-PRESSED", "KOKUGO", "KUN",
                               "LAST-DUP", "LE", "LEVEL", "LIKE", "LINE", "LINE-COUNT", "LINE-NUMBER", "LINK",
                               "LIST", "LOWVALUES", "LQ", "LS", "LT",
                               "MASK", "MATCHED", "MEND", "MESSAGE", "MOVE", "MSTART", 
                               "NE", "NEWPAGE", "NOMASK", "NOPRINT", "NOT", "NOTE", "NOVERIFY", "NQ", "NULL", 
                               "OF", "OR", "OTHERWISE",
                               paste("PA", 1:3, sep=""), "PAGE-COUNT", "PAGE-NUMBER", "PARM-REGISTER", "PATH-ID", "PATTERN",
                               "PERFORM", "POINT", "POS", "PRIMARY", "PRINT", "PROC", "PROCEDURE", "PROGRAM", "PUT",
                               "READ", "RECORD", "RECORD-COUNT", "RECORD-LENGTH", "REFRESH", "RELEASE", "RENUM", "REPEAT", 
                               "REPORT", "REPORT-INPUT", "RESHOW", "RESTART", "RETRIEVE", "RETURN-CODE", "ROLLBACK", "ROW",
                               "S", "SCREEN", "SEARCH", "SECONDARY", "SELECT", "SEQUENCE", "SIZE", "SKIP", "SOKAKU", "SORT",
                               "SQL", "STOP", "SUM", "SYSDATE", "SYSDATE-LONG", "SYSIN", "SYSIPT", "SYSLST", "SYSPRINT", "SYSSNAP", "SYSTIME",
                               "TALLY", "TERM-COLUMNS", "TERMINATION", "TERM-NAME", "TERM-ROWS", "TITLE", "TO", "TRANSFER", "TRC",
                               "UNIQUE", "UNTIL", "UPDATE", "UPPERCASE", "USER", "USERID", 
                               "VALUE", "VERIFY",
                               "W", "WHEN", "WHILE", "WORK", "WRITE",
                               "X", "XDM", "XRST")

java.reserved.words <- c( "abstract", "continue", "for", "new", "switch", "assert", "default", "goto", "package",
                          "synchronized", "boolean", "do", "if", "private", "this", "break", "double", "implements",
                          "protected", "throw", "byte", "else", "import", "public", "throws", "case", "enum", 
                          "instanceof", "return", "transient", "catch", "extends", "int", "short", "try", "char",
                          "final", "interface", "static", "void", "class", "finally", "long", "strictfp",
                          "volatile", "const", "float", "native", "super", "while")

# Turn a string into a vector of words
# for comparability across bags of words, also strip out punctuation and
# numbers, and shift all letters into lower case
# Input: character string
# Output: vector of words (charaacter strings)
strip.cobol.text <- function(txt) {  
  
  #txt <- tolower(txt)
  
  #change other non-alphanumeric characters to spaces
  txt <- gsub("[^a-zA-Z0-9]"," ",txt)
  
  # split and make one vector
  
  txt <- strsplit(txt," ")
  
  txt <- unlist(txt)
  
  # convert to lowercase
  txt <- tryTolower(txt)
  
  # remove words of one character or less
  txt <- txt[nchar(txt) > 1]
  
  # remove numeric data
  txt <- txt[!is.numeric(txt)]
  return(txt)
}


# Turn a string into a vector of words
# for comparability across bags of words, also strip out punctuation and
# numbers, and shift all letters into lower case
# Input: character string
# Output: vector of words (charaacter strings)
strip.java.text <- function(txt, lengthLowerBound=1) {  
  
  #txt <- tolower(txt)
  
  #change other non-alphanumeric characters to spaces
  txt <- gsub("[^a-zA-Z0-9]"," ",txt)
  
  txt <- gsub('((?<=[a-z])[A-Z]|[A-Z](?=[a-z]))', " \\1", txt, perl=TRUE)
  
  # split and make one vector
  
  txt <- strsplit(txt," ")
  
  txt <- unlist(txt)
  
  # convert to lowercase
  txt <- tryTolower(txt)
  
  # remove words of one character or less
  txt <- txt[nchar(txt) > lengthLowerBound]
  
  # remove numeric data
  txt <- txt[!is.numeric(txt)]
  return(txt)
}

#


# Turn a string into a vector of words
# for comparability across bags of words, also strip out punctuation and
# numbers, and shift all letters into lower case
# Input: character string
# Output: vector of words (charaacter strings)
strip.easytrieve.text <- function(txt, lengthLowerBound=1) {  
  
  #txt <- tolower(txt)
  
  #change other non-alphanumeric characters to spaces
  txt <- gsub("[^a-zA-Z0-9]"," ",txt)
  
  # split and make one vector
  
  txt <- strsplit(txt," ")
  
  txt <- unlist(txt)
  
  # convert to lowercase
  txt <- tryTolower(txt)
  
  # remove words of one character or less
  txt <- txt[nchar(txt) > lengthLowerBound]
  
  # remove numeric data
  txt <- txt[!is.numeric(txt)]
  return(txt)
}

# Stem a vector of words using a user-defined natural language stemmer
# Input: vector of words, natural language
# Output: vector of stemmed words
stem.words <- function(x,l) {
  require(RTextTools)
  return(wordStem(x, l))
  
}

# Remove the stopwords based on a user-defined natural language
# Input: vector of words, natural language
# Output: vector of reduced words
remove.stopwords <- function(x,l) {
  require(tm)
  x <- x[! x %in% stopwords(l)]
  return(x)
}

#Precondition: the reserved words are already removed
# Prepare the list of word vectors based on a user-defined natural language
# Input: a list of vector of words, natural language
# Output: a list of vector of reduced words
prepare.natural.lang.list <- function(x, l) {
  
  x <- lapply(x, remove.stopwords, l)
  return(lapply(x, stem.words, l))  
}

# Prepare the list of word vectors based on a user-defined programming language
# Input: a list of vector of words, programming language
# Output: a list of vector of reduced words
prepare.prog.lang.list <- function(x, p) {
  
  return(lapply(x, remove.reservedwords, p))  
}


# Remove the stopwords based on a user-defined programming language
# Input: vector of words, programming language
# Output: vector of reduced words
remove.reservedwords <- function(x,p) {  
  x <- x[! x %in% reservedwords(p)]
  return(x)  
}


# Returns the reservedwords based on a user-defined programming language
# Input: programming language
# Output: reserved words
reservedwords <- function(p) {
  if (p == "cobol") return (tolower(cobol.reserved.words))
  else if (p == "java") return (tolower(java.reserved.words))
  else if (p== "pli" ) return(tolower(c(pli.reserved.words, pli.keywords)))
  else if (p== "easytrieve" ) return(tolower(easytrieve.reserved.words))
  else return(c())
}


#Create BoW from the vector of words
#Input: a vector of words
#Output: a BoW vector with named columns
make.BoW <- function(x) {
  
  BoW <- c()
  
  if(length(x) == 0)
    return(BoW)
  
  for(i in 1:length(x)){
  
    if (x[i] %in% names(BoW))
      BoW[x[i]] <- BoW[x[i]] + 1      
    else 
      BoW[x[i]] <- 1
  }
  
  #print(length(BoW))
  #names(BoW) <- names(x)
  return(BoW)
  
  
}

#Create BoW list from a list of vector of words
#Input: a list of vector of words
#Output: a list of BoW vectors
make.BoW.list <- function(x) {
  
  return (lapply(x, make.BoW))
  
}


# Rescale the columns of a data frame or array by a given weight vector
# Input: arrray, weight vector
# Output: scaled array
scale.cols <- function(x,s) {
  return(t(apply(x,1,function(x){x*s})))
}

# Rescale rows of an array or data frame by a given weight vector
# Input: array, weight vector
# Output: scaled array
scale.rows <- function(x,s) {
  return(apply(x,2,function(x){x*s}))
}

# Compute inverse document frequency weights and rescale a data frame
# Input: data frame
# Calls: scale.cols
# Output: scaled data-frame
idf.weight <- function(x) {
  # IDF weighting
  doc.freq <- colSums(x>0)
  doc.freq[doc.freq == 0] <- 1
  w <- log(nrow(x)/doc.freq)
  return(scale.cols(x,w))
}

# Normalize vectors by the sum of their entries
# Input assumed to be a set of vectors in array form, one vector per row
# QUESTION FOR STUDENTS: What is the 1e-16 doing here?
# Input: matrix/data frame/array
# Calls: scale.rows()
# Output: matrix/data frame/array
div.by.sum <- function(x) {
  scale.rows(x,1/(rowSums(x)+1e-16))
}

# Normalize vectors by their Euclidean length
# Input assumed to be a set of vectors in array form, one vector per row
# QUESTION FOR STUDENTS: What is the 1e-16 doing here?
# Input: array
# Calls: scale.rows()
# Output: array
div.by.euc.length <- function(x) {
  scale.rows(x,1/sqrt(rowSums(x^2)+1e-16))
}


# Remove columns from a ragged array which only appear in one row
# Input: Ragged array (vectors with named columns)
# Output: Ragged array, with columns appearing in only one vector deleted
remove.singletons.ragged <- function(x) {
  # Collect all the column names, WITH repetition
  col.names <- c()
  for(i in 1:length(x)) {
    col.names <- c(col.names, names(x[[i]]))
  }
  # See how often each name appears
  count <- table(col.names)
  # Loop over vectors and keep only the columns which show up more than once
  for(i in 1:length(x)) {
    not.single <- (count[names(x[[i]])] > 1)
    x[[i]] <- x[[i]][not.single]
  }
  return(x)
}


# Standardize a ragged array so all vectors have the same length and ordering
# Supplies NAs for missing values
# Input: a list of vectors with named columns
# Output: a standardized list of vectors with named columns
standardize.ragged <- function(x) {
  # Keep track of all the column names from all the vectors in a single vector
  col.names <- c()
  # Get the union of column names by iterating over the vectors - using
  # setdiff() is faster than taking unique of the concatenation, the more
  # obvious approach
  for(i in 1:length(x)) {
    col.names <- c(col.names, setdiff(names(x[[i]]),col.names))
  }
  # put the column names in alphabetical order, for greater comprehensibility
  col.names <- sort(col.names)
  # Now loop over the vectors again, putting them in order and filling them out
  # Note: x[[y]] returns NA if y is not the name of a column in x
  for (i in 1:length(x)) {
    x[[i]] <- x[[i]][col.names]
    # Make sure the names are right
    names(x[[i]]) <- col.names
  }
  return(x)
}

# Turn a list of bag-of-words vectors into a data frame, one row per bag
# Input: list of BoW vectors (x),
#   list of row names (row.names, optional),
#   flag for whether singletons should be removed,
#   flag for whether words missing in a document should be coded 0
# Output: data frame, columns named by the words and rows matching documents
make.BoW.frame <- function(x,row.names,remove.singletons=TRUE, remove.constants=TRUE,
                           absent.is.zero=TRUE) {
  # Should we remove one-time-only words?
  if (remove.singletons) {
    y <- remove.singletons.ragged(x)
  } else {
    y <- x
  }
  # Standardize the column names
  y <- standardize.ragged(y)
  # Transform the list into an array
  # There are probably slicker ways to do this
  z = y[[1]] # Start with the first row
  if (length(y) > 1) { # More than one row?
    for (i in 2:length(y)) {
      z = rbind(z,y[[i]],deparse.level=0) # then stack them
    }
  }
  # Make the data frame
  # use row names if provided
  
  print(dim(z))
  
  if(missing(row.names)) {
    BoW.frame <- data.frame(z)
  } else {
    BoW.frame <- data.frame(z,row.names=row.names)
  }
  
  print(dim(BoW.frame)[2])
  
  
  if (absent.is.zero) {
    # The standardize.ragged function maps missing words to "NA"; replace
    # those with zeroes to simplify calculation
    BoW.frame <- apply(BoW.frame,2,function(q){ifelse(is.na(q),0,q)})
  }
  
  if (remove.constants) {
    for (i in dim(BoW.frame)[2]:1){
      if(sd(BoW.frame[,i])==0)
        BoW.frame<-BoW.frame[,-i] # Get rid of columns with a standard deviation of 0  
    }
  }
  
  
  
  return(BoW.frame)
}

# Produce a distance matrix from a data frame
# Assumes rows in the data frame (or other array) are vectors
# By default uses Euclidean distance but could use other functions as well
# cf. the built-in function dist()
# Input: array, optional distance function
# Calls: sq.Euc.dist()
# Output: matrix of distances
distances <- function(x,fun) {
  # Use Euclidean distance by default
  if (missing(fun)) {
    return(sqrt(sq.Euc.dist(x)))
  }
  # otherwise, run the function fun over all combinations of rows
  else {
    # make a new array
    n <- nrow(x)
    d <- array(NA,c(n,n),list(rownames(x),rownames(x))) #preserve row-names,
    # but also make them column names
    # iterate over row-pair combinations
    for(i in 1:n) {
      for(j in 1:n) {
        # fill the entries of the array
        d[i,j] <- fun(x[i,],x[j,])
      }
    }
    # we're done
    return(d)
  }
}


# Produce a similarity matrix from a data frame
# Assumes rows in the data frame (or other array) are vectors
# By default uses cosine similarity but could use other functions as well
# cf. the built-in function dist()
# Input: array, optional similarity function
# Calls: cos.sim()
# Output: matrix of similarities
similarities <- function(x,fun) {
  # Use Euclidean distance by default
  if (missing(fun)) {
    return(cos.sim(x))
  }
 
}


# calculate the cosine similariry between two sets of vectors
# Input: vectors in matrix form (one per row),
#        second set of vectors ditto (if missing assumed equal to first)
# Output: matrix of similarities
cos.sim <- function(x,y=NULL) {
  require(lsa)
  return(cosine(x, y))
  
}

# calculate the squared Euclidean distances between two sets of vectors
# specifically, d[i,j] is the squared distance from x[i,] to y[j,]
# Input: vectors in matrix form (one per row),
#        second set of vectors ditto (if missing assumed equal to first)
# Output: matrix of distances
sq.Euc.dist <- function(x,y=x) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  nr=nrow(x)
  nc=nrow(y)
  x2 <- rowSums(x^2)
  xsq = matrix(x2,nrow=nr,ncol=nc)
  y2 <- rowSums(y^2)
  ysq = matrix(y2,nrow=nr,ncol=nc,byrow=TRUE)
  xy = x %*% t(y)
  d = xsq + ysq - 2*xy
  if(identical(x,y)) diag(d) = 0
  d[which(d < 0)] = 0
  return(d)
}

# For each vector (row) in matrix A, return the index of and the distance to
# the index of the closest point to the matrix B
# If the matrix B is omitted, then it's assumed to be A, but no point is
# allowed to be its own own closest match
# A pre-computed distance matrix is an optional argument, otherwise it's
# computed in the squared Euclidean metric
# Input: matrix A, matrix B (optional),
# matrix of distances between them (optional)
# Output: list of vectors, one giving the
# indices, the other giving the distances
nearest.points <- function(a,b=a,d=sqrt(sq.Euc.dist(a,b))) {
  # "allocate" a vector, giving the distances to the best matches
  b.dist = numeric(nrow(a))
  
  if (identical(a,b)) {
    diag(d) = Inf
  }
  b.which = apply(d,1,which.min)
  for (i in 1:nrow(a)) {
    b.dist[i] = d[i,b.which[i]]
  }
  return(list(which=b.which,dist=b.dist))
}

#precondition: rownames(mydata) == colnames(mydata)
normalize.names <- function(mydata) {
  names <- rownames(mydata)
  
  names <- lapply(names, function(s) {
    start = which(strsplit(s, '')[[1]]=='/')
    end = which(strsplit(s, '')[[1]]=='.')
    substring(s, start+1, end-1)
  })
  
  names <- lapply(names, function(s) {
    #if (nchar(s) > 5) return (substr(s, 1, 5))
    #else return(s)
    substr(s, 1, 5)
  })
  
  rownames(mydata) = colnames(mydata) = names
  
  return(mydata)
}


#Normalizing to range 0 to 1
#Usage: apply(m, 2, range01)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


compute_bow_identifier_name_similarity <- function(identifierNames,  nat.langs="english", prog.langs="java"){
  
  mydata <- lapply(identifierNames, function(identifierName) strip.java.text(identifierName))
  
  for (i in 1:length(prog.langs))
    mydata <- prepare.prog.lang.list(mydata, prog.langs[i])
  
  for (i in 1:length(nat.langs)) 
    mydata <- prepare.natural.lang.list(mydata, nat.langs[i])
  
  names(mydata) <- identifierNames
  
  mydata <- mydata[-which(unlist(lapply(mydata, function(data) length(data) ==0)))]
  
  mydata.BoW.list <- make.BoW.list(mydata)
  
  mydata.BoW.frame <- make.BoW.frame(mydata.BoW.list, names(mydata.BoW.list))
  
  mydata.BoW.idf.frame <- idf.weight(mydata.BoW.frame)
  
  sim <- compute_cosine_kernel(mydata.BoW.idf.frame)
  
  return(sim)
}


apply.bow <- function(dirname, pattern, nat.langs, prog.langs) {
  
  mydata <- read.directory(dirname, pattern)
  
  ################
  # Preprocessing data #
  ################
  
  for (i in 1:length(prog.langs))
    mydata <- prepare.prog.lang.list(mydata, prog.langs[i])
  
  for (i in 1:length(nat.langs)) 
    mydata <- prepare.natural.lang.list(mydata, nat.langs[i])
   
  mydata.BoW.list <- make.BoW.list(mydata)
  
  mydata.BoW.frame <- make.BoW.frame(mydata.BoW.list, names(mydata.BoW.list))
  
  mydata.BoW.idf.frame.t <- t(na.omit(mydata.BoW.idf.frame))
  
  sim <- compute_cosine_kernel(Phi_d)
  
  return(mydata.BoW.frame)
}


apply.lsa <- function(dirname, pattern, nat.langs, prog.langs) {
  ################
  # Loading data #
  ################
  
  mydata <- read.directory(dirname, pattern)
  
  ################
  # Preprocessing data #
  ################
  
  for (i in 1:length(prog.langs))
    mydata <- prepare.prog.lang.list(mydata, prog.langs[i])
  
  for (i in 1:length(nat.langs)) 
    mydata <- prepare.natural.lang.list(mydata, nat.langs[i])
  

  
  mydata.BoW.list <- make.BoW.list(mydata)
  
  mydata.BoW.frame <- make.BoW.frame(mydata.BoW.list, names(mydata.BoW.list))
  
  print(dim(mydata.BoW.frame))
  
  write.table(mydata.BoW.frame, file ="mydata-BoW-matrix.csv",row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  mydata.BoW.idf.frame <- idf.weight(mydata.BoW.frame)
  
  ################
  # Compute Dissimilarity Matrix #
  ################
  
  #Normalize by euclidean distance
  
  mydata.Euc.dist.matrix <- div.by.euc.length(mydata.BoW.idf.frame)
  
  
  #compute squared euclidean distance
  
  mydata.Euc.dist.matrix <- distances(mydata.Euc.dist.matrix)
  
  
  #compute cosine similariry
  
  mydata.BoW.idf.frame.t <- t(na.omit(mydata.BoW.idf.frame))

  
  mydata.Cos.sim.matrix <- similarities(mydata.BoW.idf.frame.t)
  
  ################
  # Writing data #
  ################
  
  
  write.table(mydata.BoW.idf.frame, file ="mydata-idf-BoW-matrix.csv",row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  write.table(mydata.Euc.dist.matrix, file ="mydata-Euc-distance-matrix.csv",row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  write.table(mydata.Cos.sim.matrix, file ="mydata-Cos-similarity-matrix.csv",row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
}


################
# Visualizing data #
################
visualize.heatmap <- function(mydata.matrix, clust.method="complete", filename) {
  require(RColorBrewer)
  require(gplots)
  
  #cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(256)
  
  distCor <- function(x) as.dist(1-cor(t(x)))
  hclustComplete <- function(x) hclust(x, method=clust.method)
  
  png("cosineSimHeatmap.png", pointsize = 15, width = 3840, height = 3840)
  
  #h <- heatmap.2(mydata.Euc.dist.matrix, trace="none", scale="row", zlim=c(-3,3), reorder=TRUE,
  #          distfun=distCor, hclustfun=hclustComplete, col=rev(cols), symbreak=FALSE)
  
  h <-heatmap.2(mydata.matrix, trace="none", scale="row", zlim=c(-3,3), reorder=TRUE,
                 distfun=distCor, hclustfun=hclustComplete, col=rev(cols), symbreak=FALSE)

  #dev.off()
}

extract.information <- function(prname) {
  
  
  setwd("~/workspace")
  
  #Load and compute the lexical sim matrix
  lexsim <- read.table(paste("benchmark", prname , "mydata-Cos-similarity-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  


  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
  #cfg <- read.table("benchmark/jedit-5.1.0/cfg.csv", sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  
  #Build mydata
  mydata = list(cfg=as.matrix(cfg), lexsim=as.matrix(lexsim))
  mydata <- make.compatible(mydata)
  
  LOC = 0
  
  files = rownames(mydata$cfg)
  
  for (i in 1:length(files)) {
    lines = readLines(paste("benchmark", prname ,files[i], sep="/"))
    #lines <- lines[which((lines != "") && (lines != "\t") && (lines != " "))]
    lines <- lines[which(lines != "")]
    #print(lines)
    #print(length(lines))
    
    LOC <- LOC + length(lines)
  }
  
  
  print("Number of LOC")
  print(LOC)
  
  print("Number of all classes:")
  print(dim(lexsim)[1])
  
  print("Number of included classes:")
  print(dim(mydata$cfg)[1])
  
}

unify <- function(bow1, bow2){
  
#   if(dim(bow1)[1] != dim(bow2)[1])
#     stop("different number of rows")
  
  all.names <- union(colnames(bow1), colnames(bow2))
  
#   union.bow <- lapply(all.names, function(n) sum(bow1[,n], bow2[,n]))
  
  union.bow <- matrix(0, nrow=dim(bow1)[1] + dim(bow2)[1], ncol = length(all.names), 
                      dimnames = list(c(rownames(bow1), rownames(bow2)), all.names))

  for(dname in rownames(bow1))
    for(term in colnames(bow1))
      union.bow[dname, term] <- bow1[dname, term]
    
print("done1")

for(dname in rownames(bow2))
  for(term in colnames(bow2))
    union.bow[dname, term] <- bow2[dname, term]

  return(union.bow)
}
