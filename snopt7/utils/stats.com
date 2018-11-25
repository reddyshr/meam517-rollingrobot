#! /bin/csh -f

 set PROBSDIR=$HOME/snopt7/utils
 set FILE = $1

\rm	fort.4 >& /dev/null

egrep '(Begin|BEGIN|EXIT --|No. of it|No. of maj|No. of call|constraint violn|e for solv|  Linear |Problem name)'  {$FILE}.out  >  fort.4

 $PROBSDIR/stats

 mv fort.4	{$FILE}.inn 	>& /dev/null
 mv fort.7	{$FILE}.tex	>& /dev/null
 mv fort.9	{$FILE}.sum	>& /dev/null

#\rm {$FILE}.inn >& /dev/null
#\rm {$FILE}.out >& /dev/null
