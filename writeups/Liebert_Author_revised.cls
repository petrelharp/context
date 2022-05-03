\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{Liebert_Author_R}[2020/22/09 v1.00
^^J LaTeX class file for Liebert Journals applicable for preparing manuscript at author's end]

\LoadClass[twoside]{article}

\usepackage{mathtools,amsthm,amsfonts,amssymb}%
\usepackage[%
    natbib=true,%
    citestyle=authoryear-comp,%
    bibstyle=authortitle,%
    backend=biber,%
    date=year,%
    url=false,%
    eprint=false,%
    isbn=false,%
    giveninits=true,%
    terseinits=true,%
    uniquename=false,%
    uniquelist=false,%
    maxnames=3,%
    minnames=3%
]{biblatex}%
\usepackage[doublespacing]{setspace}%
\usepackage{threeparttable,booktabs}%
\usepackage{graphicx}%
\usepackage[figuresright]{rotating}%
\usepackage{enumitem}%

\newcommand{\keywords}[1]{\bgroup\noindent\textbf{Keywords:}
#1\egroup}%

\renewenvironment{abstract}{%
\begin{quote}\textbf{Abstract:} }%
{\end{quote}}%

\usepackage[hidelinks]{hyperref}%

\usepackage[ruled]{algorithm2e}%

\usepackage{mathrsfs}%

\newtheorem{theorem}{Theorem}%
\newtheorem{lemma}{Lemma}%
\newtheorem{corollary}{Corollary}%
\newtheorem{remark}{Remark}%

\newtheorem*{theorem*}{Theorem}%
\newtheorem*{lemma*}{Lemma}%
\newtheorem*{corollary*}{Corollary}%
\newtheorem*{remark*}{Remark}%

\AtBeginDocument{\pagestyle{headings}}%

\newcommand{\notetocomp}[1]{\bgroup\parindent\z@\raggedright%
\marginpar{\footnotesize#1\par}
\egroup}%

\renewbibmacro*{issue+date}{%
    \usebibmacro{date}%
  \newunit%
}%

\renewbibmacro*{journal+issuetitle}{%
  \usebibmacro{journal}%
  \setunit*{\addspace}%
  \iffieldundef{series}
    {}
    {\newunit
     \printfield{series}%
     \setunit{\addspace}}%
  \setunit{\addspace}%
  \usebibmacro{issue+date}%
  \setunit{\addsemicolon}%
  \printfield{volume}%
  \printfield{number}%
  \setunit{\addcolon}%
  \newunit%
}%

\renewbibmacro*{note+pages}{%
  \printfield{note}%
  \setunit{\addcolon}%
  \printfield{pages}%
  \newunit%
}%

\renewbibmacro*{doi+eprint+url}{%
  \setunit{\addsemicolon\addspace}%
    \printfield{doi}%
}%

\DeclareFieldFormat*{title}{#1}%
\DeclareFieldFormat{pages}{#1}%
\DeclareFieldFormat{doi}{doi:#1}%
\DeclareFieldFormat[article]{number}{\mkbibparens{#1}}
\DeclareNameAlias{sortname}{last-first}%
\endinput
