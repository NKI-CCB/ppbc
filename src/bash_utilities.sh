#Return file name (and the path to the file) for every file of type .php|.html|.js 
#containing the case-insensitive string "document.cookie" | "setcookie"
#r = recursive
#i = case-insensitive

egrep -ir --include=*.{php,html,js} "(document.cookie|setcookie)"