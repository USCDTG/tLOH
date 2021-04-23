\name{removeOutlierFromCalc}
\alias{removeOutlierFromCalc}
\alias{dataframe}
\alias{cols}
\alias{rows}
\alias{newValue}

\title{
Removes outliers
}
\description{
Take rows with a total count greater than 2000 and sets to NA
}

\arguments{
  \item{dataframe}{
input dataframe
}
  \item{cols}{
which column
}
  \item{rows}{
which row
}
  \item{newValue}{
what to replace
}
}

\value{
Dataframe returned
}

\author{
Michelle Webb
}

\examples{
test <- data.frame(TOTAL=c(2000,20,20,20,20,20))
removeOutlierFromCalc(test,"TOTAL",test[test$TOTAL > 2000,],NA)
}
