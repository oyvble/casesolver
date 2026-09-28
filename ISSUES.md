
Known issues (indicated by GPT-6 Astra Medium):
- Editing or deleting references through the reference editor leaves existing results active.
	-> In gui.R, f_addref, lines 2196–2250, changes are saved to refDataTABLE without updating or invalidating comparisons, fits, or match status. Your deletion fix covers “Selected profiles → Delete”, but the editor provides another deletion route. Changing alleles also leaves LRs calculated from the previous genotype visible alongside the new genotype.

- Import errors are silently discarded, while the console reports success.
	-> n gui.R, lines 1007–1008, error = function(e) e returns the error object, but nothing captures or prints it. Execution then prints "Imported from file: ...". Consequently, a failed import looks successful. Errors occurring after some assignments can also leave partially imported content.

Suggestive updates: 

 - Modify unknown names in the report (manually select name).
 
 - Use the newer calcMLE function directly instead of the depricated contLikMLE function from euroformix R-package.
 
 - getStructuredData: Handle ambiguous grouping of partial single-source profiles. Even with minLoc = 10, A–C and B–C can match while A–B disagrees. The current approach can overwrite unknown assignments and produce order-dependent consensus profiles.

 - Issue with WoE module: resList object from get("resWOEeval",envir=nnTK) expects fitted model as indexes.
   However, the insertion of table resList$resTable together with indices may be problematic!

 - Don't re-caclulate all comparisons after clicking "compare" when including new references. 
	- Only calculate new 'EVID~POI|(COND,NOC)' combinations
	- Qualitative, Quantitative 
  
 - Store name of deselected profile in report.
 
 - Profile manipulation:
	- Collapse similar Refs (from IBS). Useful for reducing number of unknowns.
	
 - Possible to change alligning in table values and headers (to left aligning)? 

 - When using EFMex as an option in WoE: Conditional(s) expanded with fitted?  

 - Potential Bugs (issues):
	- When doing DC after single quanLR calculations and REFERENCE is removed.
	- Identical references not removed if different names? Example "Case 7154"
	- When creating report with QuanLR matchlist but QuanLR not selected as model.  
	- When trying to create a report, getting this error in R: Error in plot.window(…) : need finite ‘xlim’ values.
  
 - Known (fixed) issues: 
	- Crashing when DC with 1 contribution, condition on known contributor with missing markers (caused by EFM v3.0.4 or earlier).

CRASH WHEN
- PlotTopEPG/MPS
- Selecting too many to condition on (DC). MatchList panel.
- ISSUE: CANT SHOW MPSplots in report (when LUS+)

ISSUE: An error occured after running WoE evidence calculations...
Error in structure(.External(.C_dotTclObjv, objv), class = "tclObj") :
  [tcl] bad window path name ".8.5.1.1.1.3".
