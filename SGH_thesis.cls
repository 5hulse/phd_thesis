% SGH_thesis
% Written by Simon Hulse (simon.hulse@chem.ox.ac.uk), 2021

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{SGH_thesis}[2021/01/14 Custom thesis class by Simon Hulse]

\LoadClass[openright,twoside,12pt,a4paper]{report}

\newif\ifelectronicversion
\electronicversiontrue

\DeclareOption{bindversion}{\electronicversionfalse}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions\relax

%%% colours
\RequirePackage{xcolor}
\definecolor{boxcolor}{HTML}{000000}


%%% font configuration
\RequirePackage{fontspec}
\RequirePackage[%
    math-style=ISO,%
    bold-style=ISO,%
    nabla=upright,%
]{unicode-math}

\setmainfont{EBGaramond}[%
    Path           = fonts/EBGaramond/,%
    Extension      = .otf,%
    UprightFont    = *-Regular,%
    BoldFont       = *-SemiBold,%
    ItalicFont     = *-Italic,%
    BoldItalicFont = *-BoldItalic,%
    FontFace       = {sb}{n}{*-SemiBold},%
]
\setmathfont{Garamond-Math}[%
    Path           = fonts,%
    Extension    = .otf,%
    Path         = fonts/,%
    StylisticSet = {2,4,8},%
]
\setmonofont{FiraMono}[%
    Scale       = MatchLowercase,%
    Path        = fonts/FiraMono/,%
    Extension   = .ttf,%
    UprightFont = *-Regular,%
    BoldFont    = *-Bold,%
]


\RequirePackage{eso-pic}
\RequirePackage{graphicx}

% 1.5 line spacing
\linespread{1.25}

% college name
\def\college#1{\gdef\@college{#1}}
% degree name
\def\degree#1{\gdef\@degree{#1}}
% year and term
\def\degreedate#1{\gdef\@degreedate{#1}}
% supervisors
\def\supervisors#1{\gdef\@supervisors{#1}}
\newif\ifsupervisors

\ifelectronicversion
  \RequirePackage[%
  includehead,%
  top      = 2.5cm,%
  bottom   = 2.5cm,%
  inner    = 2.8cm,%
  outer    = 2.8cm,%
  headsep  = 0.8cm,%
  footskip = 1.5cm,%
  ]{geometry}
\else
  \RequirePackage[%
  includehead,%
  top      = 2.5cm,%
  bottom   = 2.5cm,%
  inner    = 3.3cm,%
  outer    = 2.3cm,%
  headsep  = 0.8cm,%
  footskip = 1.5cm,%
  ]{geometry}
\fi

%%% Headers
\RequirePackage{fancyhdr}
\pagestyle{fancy}

\renewcommand{\chaptermark}[1]{\markboth{{\bfseries \thechapter}\hspace{10pt}#1}{}}
\renewcommand{\sectionmark}[1]{\markright{{\bfseries \thesection}\hspace{10pt}#1}}

% Headers and footers for "plain" pages: essentially chapter title pages.
\fancypagestyle{plain}{%
    \fancyhf{}%
    \fancyhead[LE]{%
    \makebox[0pt][r]{%
    \colorbox{boxcolor}{\makebox[\textwidth][r]{\textcolor{white}{\thepage}\enskip}}\hspace*{1em}}%
  }
    \fancyhead[RO]{%
    \makebox[0pt][l]{%
      \hspace*{1em}\colorbox{boxcolor}{\makebox[\textwidth][l]{\enskip\textcolor{white}{\thepage}}}}%
    }
}%

% Headers and footers for "typical" pages.
\fancypagestyle{fancy}{%
    \fancyhf{}%
    \fancyhead[LE]{%
    \makebox[0pt][r]{%
      \colorbox{boxcolor}{\makebox[\textwidth][r]{\fontseries{sb}\selectfont\textcolor{white}{\thepage}\enskip}}\hspace*{1em}}%
    \headrulefill\hspace*{1em}{\leftmark}\hspace*{1em}\headrulefill%
    }%
    \fancyhead[RO]{%
    \mbox{}\headrulefill \hspace*{1em}{\rightmark}\hspace*{1em}\headrulefill%
    \makebox[0pt][l]{%
      \hspace*{1em}\colorbox{boxcolor}{\makebox[\textwidth][l]{\enskip \fontseries{sb}\selectfont\textcolor{white}{\thepage}}}}%
    }
}

\renewcommand\headrulewidth{0pt}
\def\headrulefill{\leaders\hrule width 0pt height 3pt depth -2.8pt \hfill}

% clear empty trailing left page
% http://www.texfaq.org/FAQ-reallyblank
\let\origdoublepage\cleardoublepage
\newcommand{\clearemptydoublepage}{%
  \clearpage
  {\pagestyle{empty}\origdoublepage}%
}
\let\cleardoublepage\clearemptydoublepage

%%% Paragraph skip
\RequirePackage{parskip}
%%% Chapter and Section titles
\RequirePackage[explicit]{titlesec}
\RequirePackage{anyfontsize}
\RequirePackage{setspace}

\titleformat{\chapter}[hang]{}{%
  \thispagestyle{empty}
  \parbox{0.8\textwidth}{
      \raggedright\bfseries\fontsize{28}{36}\selectfont#1
  }
  \hspace*{3em}
  \colorbox{boxcolor}{%
    \begin{minipage}[c][4em][c]{.5\textwidth}
      {\hspace*{1em}\fontsize{50}{48}\selectfont\bfseries\textcolor{white}{\thechapter}}
    \end{minipage}
}}{0pt}{}[\vspace*{4em}]

\renewcommand*\contentsname{Table of Contents}
\titleformat{name=\chapter,numberless}[display]
    {}
    {}
    {-80pt}
    {\centering{\bfseries\fontsize{28}{36}\selectfont #1}}

%%% ToC formatting
\RequirePackage{titletoc}
% \tocless will be used to make summary not present in ToC
\newcommand{\nocontentsline}[3]{}
\newcommand{\tocless}[2]{\bgroup\let\addcontentsline=\nocontentsline#1{#2}\egroup}

\titlecontents{chapter}
[0pt]
{\vspace{15pt}\large}
{\bfseries\contentslabel{0pt}\hspace{20pt}}
{\bfseries\contentslabel{0pt}}
{\hfill\upshape\bfseries\contentspage}
[]

\titlecontents{section}
[20pt]
{\vspace{5pt}\bfseries}
{\contentslabel{0pt}\hspace{30pt}}
{}
{\mdseries\titlerule*[.5pc]{.}\bfseries\contentspage}
[]

\titlecontents{subsection}
[40pt]
{}
{\contentslabel{0pt}\hspace{35pt}}
{\contentsmargin{0pt}\bfseries}
{\titlerule*[.5pc]{.}\contentspage}
[]

\titlecontents{subsubsection}
[60pt]
{}
{}
{\contentsmargin{0pt}\itshape}
{\titlerule*[.5pc]{.}\contentspage}
[]

\titleformat{\paragraph}[runin]{\itshape}{}{0pt}{#1}[]
\titlespacing{\paragraph}{0pt}{10pt}{10pt}[]

%% Caption formatting
\RequirePackage[%
    font=singlespacing,%
    labelfont=small,bf,%
    textfont=small,%
]{caption}

\renewcommand*{\maketitle}{%
\begin{titlepage}
  \AddToShipoutPictureBG*{%
  	\AtTextCenter{\makebox(0,0){\includegraphics[width=\textwidth]{theme_figures/beltcrest_grey.pdf}}%
  }}
	\begin{center}
  	\begin{minipage}[t][17.2cm][b]{\textwidth}
    	\begin{center}
        \huge
        \@title\\[1cm]
      	\bfseries\LARGE{\@author}\\[0.5cm]
      	\large\mdseries\upshape
      	\@college\\
      	University of Oxford\\[1cm]
        A thesis submitted for the degree of\\
        \@degree\\
        \@degreedate\\[1cm]
        \ifsupervisors
      	  Supervisors:\\
          \itshape
          \@supervisors
        \fi
    	\end{center}
  	\end{minipage}
	\end{center}
\end{titlepage}
}

\makeatletter
\newenvironment{chapterquote}[1]
    {%
        \setlength{\@tempdima}{2em}%
        \def\chapquote@author{#1}%
        \vspace{-40pt}
        \parshape 1 \@tempdima \dimexpr\textwidth-2\@tempdima\relax\itshape
    }
    {%
        \par\normalfont\hfill--\ \chapquote@author\hspace*{\@tempdima}\par\bigskip\bigskip
    }
\makeatother

\newenvironment{dedication}{%
    \cleardoublepage
    \thispagestyle{empty}
    \begin{center}
    \vspace*{1.5cm}
    \end{center}
    \vspace{0.5cm}
    \begin{verse}\begin{center}}%
    {\end{center}\end{verse}}

\newenvironment{frontmattersection}[1]{%
    \cleardoublepage
    \newpage
    \begin{center}
    {\bfseries\fontsize{28}{36}\selectfont #1\\\vspace{1.5em}}
    \end{center}}%
    {}

\renewenvironment{abstract}
{\begin{frontmattersection}{Abstract}}
{\end{frontmattersection}}

\newenvironment{acknowledgements}
{\begin{frontmattersection}{Acknowldegements}}
{\end{frontmattersection}}

\newenvironment{acronyms}
{\begin{frontmattersection}{Acronyms}}
{\end{frontmattersection}}

\RequirePackage{enumitem}
\setlist[itemize,1]{%
    label={\raisebox{.45ex}{\tiny$\blacksquare$}},%
    labelwidth=0pt,%
    leftmargin=2em,%
    itemsep=-5pt,%
}
\setlist[enumerate,1]{
    label={\arabic*.},%
    leftmargin=2em,%
    itemsep=-5pt,%
}

%%% Bibliography configuration
\RequirePackage[
    doi=false,
    isbn=false,
    url=false,
    eprint=false,
    sorting=none,
    style=numeric-comp,
    giveninits=true,
    maxbibnames=20,
]{biblatex}
\DeclareDelimFormat{postnotedelim}{\addcolon\space}
\DeclareFieldFormat{postnote}{\mknormrange{#1}}

%%% Footnotes
% Symbols instead of numbers
\RequirePackage[symbol]{footmisc}
\renewcommand{\thefootnote}{\fnsymbol{footnote}}
% No footnote indent
\makeatletter
\renewcommand\@makefntext[1]{%
  \noindent{\@makefnmark}#1}
\makeatother

%%% Hyperlink configuration
\usepackage{hyperref}
\definecolor{linkblue}{HTML}{0000ee}
\hypersetup{
    colorlinks=true,
    linkcolor=black,
    urlcolor=linkblue,
    citecolor=black,
}
\def\UrlBreaks{\do\/\do-\do.}

%%% Algorithm configuration
\RequirePackage{algorithm,algpseudocode}
\numberwithin{algorithm}{chapter}
\makeatletter
\algrenewcommand\ALG@beginalgorithmic{\footnotesize}
\makeatother
\captionsetup[algorithm]{%
    font=singlespacing,%
    labelfont={small,bf},%
    textfont=small,%
}

%%% Internal reference configuration
\RequirePackage[capitalize,noabbrev]{cleveref}
\crefformat{equation}{Equation~#2#1#3}
\Crefformat{equation}{Equation~#2#1#3}
\crefrangeformat{equation}{Equations~(#3#1#4) to~(#5#2#6)}
\Crefrangeformat{equation}{Equations~(#3#1#4) to~(#5#2#6)}
\crefmultiformat{equation}{Equations~#2#1#3}%
{ and~#2#1#3}{, #2#1#3}{ and~#2#1#3}
\Crefmultiformat{equation}{Equations~#2#1#3}%
{ and~#2#1#3}{, #2#1#3}{ and~#2#1#3}

\makeatletter
\crefformat{tcb@cnt@mylisting}{Code Listing~#2#1#3}
\Crefformat{tcb@cnt@mylisting}{Code Listing~#2#1#3}
\crefrangeformat{tcb@cnt@mylisting}{Code Listings~(#3#1#4) to~(#5#2#6)}
\Crefrangeformat{tcb@cnt@mylisting}{Code Listings~(#3#1#4) to~(#5#2#6)}
\crefmultiformat{tcb@cnt@mylisting}{Code Listings~#2#1#3}%
{ and~#2#1#3}{, #2#1#3}{ and~#2#1#3}
\Crefmultiformat{tcb@cnt@mylisting}{Code Listings~#2#1#3}%
{ and~#2#1#3}{, #2#1#3}{ and~#2#1#3}
\makeatother

%%% Code listing configuration
\RequirePackage{minted,tcolorbox,listings}
\tcbuselibrary{listings,minted,breakable,skins}
\definecolor{code_bg}{HTML}{f0f0f0}
\renewcommand{\theFancyVerbLine}{\ttfamily\textbf{\small
    \textcolor{black}{\arabic{FancyVerbLine}}}}
\newtcbinputlisting[
    list inside=mylist,
    auto counter,
    number within=chapter,
]{\mylisting}[5]{%
    listing only,
    minted language={#1},
    listing file={#2},
    list entry={\protect\numberline{\thetcbcounter} #3},
    title={{\small\textbf{Code Listing \thetcbcounter:} #4}},
    label={#5},
    boxrule=0pt,
    titlerule=0.5pt,
    borderline north={0.5pt}{0pt}{black},
    borderline south={0.5pt}{-0.5pt}{black},
    code={\setlength{\parskip}{0.05in}\singlespacing},
    colframe=black,
    colback=code_bg,
    coltitle=black,
    colbacktitle=white,
    arc=0pt,
    outer arc=0pt,
    enhanced,
    breakable,
    top=0pt,
    bottom=0pt,
    width=\textwidth,
    minted options = {
        style=friendly,
        linenos,
        tabsize=0,
        breaklines,
        autogobble=false,
        mathescape,
        escapeinside=||,
        numbersep=2.2mm,
        xleftmargin=3mm
    },
    overlay = {
    \begin{tcbclipinterior}
        \fill[code_bg] (frame.south west) rectangle ([xshift=5mm]frame.north west);
    \end{tcbclipinterior}
    }
}
\renewcommand{\contentsname}{Contents}
\renewcommand{\listfigurename}{Figures}
\renewcommand{\listtablename}{Tables}
\renewcommand{\listalgorithmname}{Algorithms}
