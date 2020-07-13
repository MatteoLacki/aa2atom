from aa2atom import aa2atom
import re

aaseq = "AAPAPAAPGGL"
aa2atom(aaseq)
aa2atom(aaseq, no_water=)

A1 = aa2atom("AAPAPAAPGGL")
A2 = aa2atom("AASSAPAAPGGL")

A1 + A2

W = Counter("AAPAPAAPGGL")

# check these calculations
# get rid of the LC dependency

