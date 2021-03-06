% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/textEditor.R
\name{textEditor}
\alias{textEditor}
\title{textEditor}
\usage{
textEditor(env)
}
\arguments{
\item{env}{an environment with elements items1 and items2}
}
\description{
User select items between two tables
}
\details{
The function returns the items selected in the two lists (when window is closed)
}
\examples{
\dontrun{
env = new.env( parent = emptyenv() ) 
assign("text","this is my text",envir=env)
textedior(env)
}
}
