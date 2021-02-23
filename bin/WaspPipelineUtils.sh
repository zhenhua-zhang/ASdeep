#!/bin/bash
# Utilities

# ENVs {
export ASEPL_VERSION
ASEPL_VERSION="0.1.0"
ASEPL_FILENAME=$(basename "$0")
# } // ENVs

# Text color and format {
if [[ -x /usr/bin/tput ]] && tput setaf 1 &> /dev/null; then
    tput sgr0
    fmtBold=$(tput bold)
    fmtReset=$(tput sgr0)
    fmtItalic=$(tput sitm)

    fmtRed=$(tput setaf 1)
    fmtBlue=$(tput setaf 4)
    fmtCyan=$(tput setaf 6)
    fmtGrey=$(tput setaf 8)
    fmtBlack=$(tput setaf 0)
    fmtGreen=$(tput setaf 2)
    fmtWhite=$(tput setaf 7)
    fmtYellow=$(tput setaf 3)
    fmtPurple=$(tput setaf 5)
else
    fmtBold="\e[1m"
    fmtReset="\e[0m"
    fmtItalic="\e[3m"

    fmtRed="\e[1;31m"
    fmtBlue="\e[1;34m"
    fmtCyan="\e[1;36m"
    fmtGrey="\e[1;37m"
    fmtBlack="\e[1;30m"
    fmtGreen="\e[1;32m"
    fmtWhite="\e[1;37m"
    fmtYellow="\e[1;33m"
    fmtPurple="\e[1;35m"
fi
# } // Text color and format


# Error information, exit -1
echoErro() {
    echo -e "$fmtBold$fmtRed[E]: $1$fmtReset" >&2 && exit -1
}

# Warning information
echoWarn() {
    echo -e "$fmtBold$fmtYellow[W]:$fmtReset $*" >&2
}

# General information
echoInfo() {
    echo -e "$fmtBold$fmtWhite[I]:$fmtReset $*"
}

echoVersion() {
    cat << EOF

$ASEPL_FILENAME, Version ${ASEPL_VERSION:=UNKNOWN}

EOF
}

# Echo help for the script
echoHelp() {
    cat <<EOF

$ASEPL_FILENAME, Version ${ASEPL_VERSION:=UNKNOWN}

It is recommended to use the a configuration file instead of commandline
arguments and flags. Please check the \`./config\` file for more information.

Help:
  -c, --conf     Required.
      The confituration file to supply values for the variables
  -h, --help     Optional.
    Print this help context and exit.
  -v, --version  Optional. Action: store_true
    Print version of current script and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
}
