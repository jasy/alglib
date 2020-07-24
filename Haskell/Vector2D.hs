#!/usr/bin/env runghc
dot(x1,y1)(x2,y2)=x1*x2+y1*y2
cross(x1,y1)(x2,y2)=x1*y2-x2*y1
upperArg(x,y)=y>0 || y==0 && x>0
cmpArg l r=compare(upperArg r,0)(upperArg l,cross l r)
