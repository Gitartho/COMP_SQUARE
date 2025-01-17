\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{hcmut-report}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[a4paper]{article}

\RequirePackage[utf8]{vntex,inputenc}
\RequirePackage{anyfontsize,amsmath,amssymb}
\RequirePackage[T5]{fontenc}

\RequirePackage{geometry,graphicx,indentfirst,setspace}
\RequirePackage{booktabs}

\RequirePackage{fancyhdr,lastpage,hyperref}
\RequirePackage[super,square,sort]{natbib}
\RequirePackage[nameinlink]{cleveref} % use \Cref{} instead of \ref{}

\RequirePackage[english]{babel}

\newcommand{\Month}{
  \ifcase \the\month
  \or January
  \or February
  \or March
  \or April
  \or May
  \or June
  \or July
  \or August
  \or September
  \or October
  \or November
  \or December
  \fi
}

\newcommand{\@deptname}{Department of Aerospace Engineering}
\newcommand{\odeptname}[1]{\renewcommand{\@deptname}{#1}}

\newcommand{\@coursename}{Course name}
\newcommand{\ocoursename}[1]{\renewcommand{\@coursename}{#1}}

\newcommand{\@reporttype}{Assignment Report}
\newcommand{\oreporttype}[1]{\renewcommand{\@reporttype}{#1}}

\renewcommand{\@title}{Report title}
\newcommand{\otitle}[1]{\renewcommand{\@title}{#1}}

\newcommand{\@advisor}{Advisor}
\newcommand{\oadvisor}[1]{\renewcommand{\@advisor}{#1}}

\newcommand{\@reportplace}{Chennai, Tamil Nadu}
\newcommand{\oreportplace}[1]{\renewcommand{\@reportplace}{#1}}

\newcommand{\@reporttime}{\Month \the\year}
\newcommand{\oreporttime}[1]{\renewcommand{\@reporttime}{#1}}

\def\reportlayout{
  \geometry{
    paper=a4paper,
    vmargin=2cm,
    hmargin=2cm,
    includeheadfoot=true,
    headheight=30pt
  }

  \onehalfspacing%

  % Allow page breaks inside align* environment
  % \allowdisplaybreaks{}

  % Makes a lot of things blue, avoid at all costs
  % \everymath{\color{blue}}

  % Table of contents section text
  % \AtBeginDocument{\renewcommand*\contentsname{Contents}}

  % References section text
  % \AtBeginDocument{\renewcommand*\refname{References}}

  % Numbering format
  % {equation}{section} means <section number>.<equation number>
  % {equation}{subsection} means <section number>.<subsection number>.<equation number>
  % \numberwithin{equation}{section}

  \AtBeginDocument{\renewcommand*\bibname{References}}

  \hypersetup{
    colorlinks=true,
    linkcolor=black,
    citecolor=red,
    urlcolor=blue
  }

  \pagestyle{fancy}

  \renewcommand{\headrulewidth}{0.3pt}
  \fancyhead{}
  \fancyhead[L]{
    \begin{tabular}{rl}
      \begin{picture}(25,15)(0,0)
        \put(0,-8){\includegraphics[width=8mm,height=8mm]{graphics/iitm.png}}
      \end{picture}
      
      \begin{tabular}{l}
        {\ttfamily Indian Institute of Technology, Madras} \\
        {\ttfamily  Department of Aerospace Engineering}
        \vspace{0.15cm}
      \end{tabular}
    \end{tabular}
  }
  \fancyhead[R]{
    \begin{tabular}{l}
      \tiny \bf \\
      \tiny \bf
    \end{tabular}
  }

  \renewcommand{\footrulewidth}{0.3pt}
  \fancyfoot{}
  \fancyfoot[L]{\scriptsize \ttfamily \@reporttype}
  \fancyfoot[R]{\scriptsize \ttfamily Page {\thepage}/\pageref{LastPage}}
}

\def\coverpage{
  \thispagestyle{empty}

  \vspace*{-3\baselineskip}
  \begin{center}
    Indian Institude of Technology, Madras

    \MakeUppercase{\@deptname}
    \\[5\baselineskip]
    \includegraphics[scale=0.25]{graphics/iitm.png}
    \\[3\baselineskip]
    \begin{tabular}{c}
      \textbf{\Large \@coursename}                                 \\
      {}                                                           \\
      \midrule                                                     \\
      {}                                                           \\
      {\Large \@reporttype}                                 \\
      {}                                                           \\
      \textbf{\huge \parbox[c]{0.8\linewidth}{\centering \@title}} \\
      {}                                                           \\
      \bottomrule
    \end{tabular}
    \vfill%
    \begin{Large}
        Gitartha pran Bharadwaj \\
       \textbf{ME20B073} \\
    \end{Large}
    
    \\[3\baselineskip]
    \Make{\@reportplace, \@reporttime}
  \end{center}

  \clearpage
}
