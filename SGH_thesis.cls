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
\definecolor{headingcolor}{HTML}{002f60}


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
    BoldFont       = *-Bold,%
    ItalicFont     = *-Italic,%
    BoldItalicFont = *-BoldItalic,%
]
\setmathfont{Garamond-Math}[%
    Path           = fonts,%
    Extension    = .otf,%
    Path         = fonts/,%
    StylisticSet = {2,8},%
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
      \colorbox{headingcolor}{\makebox[\textwidth][r]{\textcolor{white}{\thepage}\enskip}}\hspace*{1em}}%
  }
    \fancyhead[RO]{%
    \makebox[0pt][l]{%
      \hspace*{1em}\colorbox{headingcolor}{\makebox[\textwidth][l]{\enskip\textcolor{white}{\thepage}}}}%
    }
}%

% Headers and footers for "typical" pages.
\fancypagestyle{fancy}{%
    \fancyhf{}%
    \fancyhead[LE]{%
    \makebox[0pt][r]{%
      \colorbox{headingcolor}{\makebox[\textwidth][r]{\textcolor{white}{\thepage}\enskip}}\hspace*{1em}}%
    \headrulefill\hspace*{1em}{\leftmark}\hspace*{1em}\headrulefill%
    }%
    \fancyhead[RO]{%
    \mbox{}\headrulefill \hspace*{1em}{\rightmark}\hspace*{1em}\headrulefill%
    \makebox[0pt][l]{%
      \hspace*{1em}\colorbox{headingcolor}{\makebox[\textwidth][l]{\enskip\textcolor{white}{\thepage}}}}%
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

%%% Chapter and Section titles
\RequirePackage[explicit]{titlesec}
\RequirePackage{anyfontsize}
\RequirePackage{setspace}

\titleformat{\chapter}[hang]{}{%
  \thispagestyle{empty}
  % \begin{minipage}[c][3em][c]{.8\textwidth}
  \parbox{0.8\textwidth}{
      \raggedright\color{headingcolor}\scshape\fontsize{28}{36}\selectfont#1
  }
  % \end{minipage}
  \hspace*{3em}
  \colorbox{headingcolor}{%
    \begin{minipage}[c][4em][c]{.5\textwidth}
      {\hspace*{1em}\fontsize{50}{48}\selectfont\bfseries\textcolor{white}{\thechapter}}
    \end{minipage}
}}{0pt}{}[\vspace*{4em}]

\renewcommand*\contentsname{Table of Contents}
\titleformat{name=\chapter,numberless}[display]
    {}
    {}
    {-80pt}
    {\centering\color{headingcolor}{\scshape\fontsize{30}{30}\selectfont #1}}

%%% ToC formatting
\RequirePackage{titletoc}
% \tocless will be used to make summary not present in ToC
\newcommand{\nocontentsline}[3]{}
\newcommand{\tocless}[2]{\bgroup\let\addcontentsline=\nocontentsline#1{#2}\egroup}

\titlecontents{chapter}
[0pt]
{\vspace{15pt}\large\color{headingcolor}}
{\contentslabel{0pt}\scshape\hspace{20pt}}
{}
{\hfill\upshape\contentspage}
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
{\contentsmargin{0pt}\scshape}
{\titlerule*[.5pc]{.}\contentspage}
[]

\titlecontents{subsubsection}
[60pt]
{}
{}
{\contentsmargin{0pt}\itshape}
{\titlerule*[.5pc]{.}\contentspage}
[]

%% Cpation formatting
\RequirePackage[%
    font=singlespacing,%
    labelfont=small,bf,sc,%
    textfont=small,%
]{caption}

\renewcommand*{\maketitle}{%
\begin{titlepage}
  \AddToShipoutPictureBG*{%
  	\AtTextCenter{\makebox(0,0){\includegraphics[width=\textwidth]{theme_figures/beltcrest_grey.pdf}}%
  }}
	\begin{center}
  	\begin{minipage}[t][16.5cm][b]{\textwidth}
    	\begin{center}
        \huge
        \@title\\[1cm]
      	\scshape\LARGE\bfseries\textcolor{headingcolor}{\@author}\\[0.5cm]
      	\large\mdseries\upshape
      	\@college\\
      	University of Oxford\\[1cm]
        A thesis submitted for the degree of\\
        \itshape
        \@degree\\
        \upshape
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
    {\color{headingcolor}\scshape\fontsize{30}{30}\selectfont #1\\ \vspace{0.5em}}
    \end{center}}%
    {}

\newenvironment{acknowledgements}
{\begin{frontmattersection}{Acknowldegements}}
{\end{frontmattersection}}

\newenvironment{acronyms}
{\begin{frontmattersection}{Acronyms}}
{\end{frontmattersection}}



\RequirePackage{enumitem}
\setlist[itemize,1]{label=\includegraphics{theme_figures/bullet.pdf}, itemsep=0pt}
\setlist[enumerate,1]{label=\textcolor{headingcolor}{(\roman*)},itemsep=0pt}

\RequirePackage[symbol]{footmisc}
\renewcommand{\thefootnote}{\fnsymbol{footnote}}
% No footnote indent
\makeatletter
\renewcommand\@makefntext[1]{%
  \noindent{\@makefnmark}#1}
\makeatother

