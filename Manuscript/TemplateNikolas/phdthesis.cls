\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{phdthesis}[2009/05/05 Standard LaTeX class]

\LoadClassWithOptions{book}
\RequirePackage{amsmath,amsthm}
\RequirePackage{float}
\RequirePackage{titlesec}
\RequirePackage{a4wide}
\RequirePackage{array}
\RequirePackage{hyperref}
\RequirePackage{url}
\RequirePackage{paralist}
\RequirePackage{cleveref}


\hypersetup{ 
 colorlinks,% 
 citecolor=black,% 
 filecolor=black,% 
 linkcolor=black,% 
 urlcolor=black 
} 


\setcounter{secnumdepth}{3}
\renewcommand{\thesubsubsection}{\thesubsection.\alph{subsubsection}}

\titleformat{\chapter}[display]
{\normalfont\bfseries\Huge}
{
{\chaptertitlename \Huge~\thechapter\large}}
{40pt}
{
}[
\thispagestyle{empty}
]

\titleclass{\part}{page}
\titleformat{\part}[display]
{\centering\normalfont}
{\huge \sffamily \bfseries \partname~\thepart}
{10pt}
{
\thispagestyle{empty}
\Huge \sffamily \bfseries
}

% \assignpagestyle{\chapter}{empty}
% \assignpagestyle{\part}{empty}

\titleformat{\section}{\normalfont\Large\bfseries}{\thesection}{1em}{}
\titleformat{\subsection}{\normalfont\large\bfseries}{\thesubsection}{1em}{}
\titleformat{\subsubsection}[runin]{\normalfont\normalsize\bfseries}{\thesubsubsection}{1em}{}
\titleformat{\paragraph}[runin]{\normalfont\normalsize}{\theparagraph}{1em}{}
\titleformat{\subparagraph}[runin]{\normalfont\normalsize}{\thesubparagraph}{1em}{}

\titlespacing*{\chapter}{0pt}{20pt}{80pt}
\titlespacing*{\section}{0pt}{3.5ex plus 1ex minus .2ex}{2.3ex plus .2ex}
\titlespacing*{\subsection}{0pt}{3.25ex plus 1ex minus .2ex}{1.5ex plus .2ex}
\titlespacing*{\subsubsection}{0pt}{3.25ex plus 1ex minus .2ex}{1.5ex plus .2ex}
\titlespacing*{\paragraph}{0pt}{3.25ex plus 1ex minus .2ex}{1em}
\titlespacing*{\subparagraph}{\parindent}{3.25ex plus 1ex minus .2ex}{1em}

\input{ttlps.def}

\newpagestyle{main}{%
%\setmarks{\chapter}{\section}
\headrule
\sethead[\thepage][][\ifthechapter{\upshape\sffamily \chaptername\ \thechapter.\ \chaptertitle}{\upshape\sffamily \chaptertitle}]
{\ifthesection{\sffamily {\thesection.\ \sectiontitle}}{\ifthechapter{\upshape\sffamily \chaptername\ \thechapter.\ \chaptertitle}{\upshape\sffamily \chaptertitle}}}{}{\thepage}
}

% \titlecontents{part}[2em]{\addvspace{20pt}\normalfont\large\bfseries}{\contentslabel{1.6em}}{\hspace*{-1.6em}}{\hfill\contentspage}
% \titlecontents{chapter}[1.5em]{\addvspace{11.8pt}\normalfont\normalsize\bfseries}{\contentslabel{1.3em}}{\hspace*{-1.3em}}{\hfill\contentspage}
% \dottedcontents{section}[3.8em]{\normalfont\normalsize}{2.3em}{1pc}
% \dottedcontents{subsection}[7em]{\normalfont\normalsize}{3.2em}{1pc}

\theoremstyle{plain}
\newtheorem*{theorem*}{Theorem}
\newtheorem{theorem}{Theorem}[chapter]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{propositiondefinition}[theorem]{Proposition-Definition}
\newtheorem{lemma}[theorem]{Lemma}
\Crefname{lemma}{Lemma}{Lemmas}
\newtheorem{corollary}[theorem]{Corollary}
\Crefname{corollary}{Corollary}{Corollaries}
\newtheorem{invariant}{Invariant}[chapter]
\newtheorem{assumption}{Assumption}[chapter]
\Crefname{assumption}{Assumption}{Assumptions}
\newtheorem{problem}[theorem]{Problem}
\Crefname{problem}{Problem}{Problems}
\newtheorem{conjecture}[theorem]{Conjecture}
\Crefname{conjecture}{Conjecture}{Conjectures}
\newtheorem{openproblem}[theorem]{Open problem}
\Crefname{openproblem}{Open problem}{Open problems}


\newtheorem{innercustomass}{Assumption}
\newenvironment{customass}[1]
  {\renewcommand\theinnercustomass{#1}\innercustomass}
  {\endinnercustomass}


\theoremstyle{definition}
\newtheorem{definition}{Definition}[chapter]

\theoremstyle{remark}
\newtheorem{remark}{Remark}[chapter]
\newtheorem{example}[remark]{Example}

\newcommand{\etal}{\textit{et al.}}
\newcommand{\etc}{\textit{etc}}

\newcolumntype{M}[1]{>{\centering}m{#1}}
\newcolumntype{P}[1]{>{\centering}p{#1}}
\newcolumntype{B}[1]{>{\centering}b{#1}}

\setlength{\pltopsep}{0.75ex}
\setlength{\plitemsep}{0.75ex}


%<*sample05>


