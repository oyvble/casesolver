Suggestive updates: 

 - Don't re-caclulate all comparisons after clicking "compare" when including new references. 
  --> Only calculate new 'EVID~POI|(COND,NOC)' combinations
  --> Qualitative, Quantitative 
  
 - Store name of deselected profile in report.
 - Profile manipulation:
  -- Collapse similar Refs (from IBS). Useful for reducing number of unknowns.
 - Possible to change alligning in table values and headers (to left aligning)? 
 - The user can select whether to use plotly or not (can be chosen under "Advanced" somewhere)

 - Potential Bugs (issues):
  -- When doing DC after single quanLR calculations and REFERENCE is removed.
  -- Identical references not removed if different names? Example "Case 7154"
  -- When creating report with QuanLR matchlist but QuanLR not selected as model.  
  -- When trying to create a report, getting this error in R: Error in plot.window(…) : need finite ‘xlim’ values.
  
 - Known (fixed) issues: 
  -- Crashing when DC with 1 contribution, condition on known contributor with missing markers (caused by EFM v3.0.4 or earlier).

CRASH WHEN
- PlotTopEPG/MPS
- Selecting too many to condition on (DC). MatchList panel.
- ISSUE: CANT SHOW MPSplots in report (when LUS+)

#Changes version 2.0.1 (30.03.2023)
 Major changes:
  - Fixed issue for exporting docx report files when using new version of officer (from v0.6.0): Major updates in createReport.

 Minor changes:
  - Adding item number to report layout.
  - Improving the importing of references from separate text file (made more robust). Markers names was case sensitive.
  - Fix warning issue in createReport function, when "selected" to inner function insList was a vector.
  - Included licence file.

#Changes version 2.0.0 (29.11.2022)
 Major changes:
  - This version requires euroformix v4 to perform additional calculations after importing data.
	-> The calculations are typically ~4x faster than earlier versions (both qualitative/quantiative models).
	
  - In Weight-of-Evidence (WoE) module: 
	-> In result table: User can now select multiple hypothesis sets.
	-> In specification panel: User can make all hypothesis sets unconditional (new button).
	-> If the LR exceeds the defined verbalized numbers: A prefix "at least" is introduced (LR value also modified).
	-> The BayesFactor and conservative LR calculation now utilizes the calcLRmcmc implementation in EFM v4. This enables both trace plots and extended iterations.

  - When importing data: 
    -> No longer exluding folders (now considered as a "file"). Possible to separete evid/ref on different folders.
	-> Name of functions that is imported can have any name (other than importData)
	
  - In IBS calculation (Ref-Ref comparison):
	-> The first column is now the number of missmatch. The table is now sorted with this.
	-> The IBS and number of markers used in comparisons are now the last columns.

  - In profile selection: Adding a "Select MatchStatus" feature:
    -> Enables user to modify the MatchStatus for selected evidence profiles.	

  - Now exporting all settings to EuroForMix: 
   -> Selected kit, population frequency, global or marker specific values of AT, fst and drop-in model.

  Minor changes:
  - Evids/Refs in match lists are now removed: MatchMatrix,MatchList(qual,quan,final)
  - Inlcude adj.LogLik to parameter table (last row).  
  - Modifying number of decimals showing parameter estimates (more pretty).
  - Separator for storing as csv is now semi-colon instead of comma (easier to open table in excel).
  - Fixed issue with swapped Yes/No when trying closing program
  - Fixed crash when creating a report where no single source are included (intent to show EPG for these)
 
#Changes version 1.8.2 (28.07.2021)

 Important changes:  
 - Reference profiles must now be an exact genotype match to show up in MatchStatus for single source profiles. In earlier versions the criterion has been inclusion:
   For instance a reference with homozygous genotype 7/7 would match a single source profile genotype of 6/7 because it was an inclusion. Thanks to Merete Ramse for highlighting this issue.
  -- A printout is given in the R-console to inform the user if there was still an inclusion (when the genotypes did not match).
  -- getStructuredData:L150-178: Including the exact genotype match criterion.

 - Language can now be imported using a text file as backup if reading with the readxl package fails.

#Changes version 1.8.1 (23.06.2021)

 Modifications:
 - gui:L961-966: Metadata appending separate the situation of vector and matrix case to avoid wrong appending.
 - calcWOEhyps:L352: An hypotheses set is no longer suggested if already calculated but differs in conditional references.
 - createReport:L535: The HTML-title argument when defining the document is now same as the header title.
 
 Fixed bugs:
 - createMatchNetwork:L61: Avoiding program to crash when LR=inf (or less than 1)
 - createreport:L417-422: insTable helpfunction caused crash when no rownames or colnames given
 - gui:L2459: Fixed crash when trying to change order of refs when no data given in data panel

#Changes version 1.8.0 (17.02.2021)

A new module is introduced: WoE evaluations for user-specific hypotheses.
 - Obtaining LR(mle,bayesian, conservative)
 - Verbal statements.
 - Other external results: Per-marker LRs, Parameter estimates, Model validation. 
 
Important updates:
 - Evid/Ref tables are now adjustable thanks to gpanedgroup function.
 - Reports can now be exported as word formats: doc (rtf) or docx (officer).
 - Marker names in report can now be fully customized (Report->'Set marker names')
 - Now allows that only reference profiles or no data are included in the case (e.g. for reporting purpose only).
 - Unknowns and deconvolved profiles are now labelled separately under "Extracted profile(s)". Only known references are now shown under References.
 - Detailed results from deconvoluted profiles can shown with new button and added to the report (RatioToNext=Probability of top ranked genotype divided by the probability of the second ranked genotypes).
 - Introducing compatibility mode of stored projects (Warning given if project version is outdated).
 - The table of deconvolved profiles can now also be sorted wrt 'Sample name' or 'similarity' (performs hierarchical clustering).
 - Options for rare alleles: The user can now set minimum frequency of imputed alleles or select whether to normalize frequencies after imputing. Options are also stored together with selected population frequency data.
 - Tooltips added to buttons in the GUI. Hover with mouse to obtain helptext (can be edited for different languages).
 - Profile Sorting assigned by the user are now remembered throughout the session (also in report): This regards the GUI tables from panels "Data", "MatchList(Qual/Quan)", "Matches".
 - "Change language" and "Change View" does not cause a full restart (all information in project is restored).

Small changes:
 - Function calcHp changed to calcMLE, and is used for quantiative LR calculation.
 - Function tabToListRef now either returns zero (empty) or two alleles (homozygous) for situations where reference has only one allele.
 - Added in getStructured:L203: Print warning if inconsistent length between number of alleles and number of peaks
 - Only "base part" (path and extension removed) of selected population file is shown in report.

Added functions: 
 - number2word: Converts large numbers to large verbal equivalent (language dependent)
 - calcWOEhyps: Wrapper function to define and calculate LR for different hypotheses.
 - calcCONS: Helpfunction to calculate Bayesian based LR and Conservative LR using MCMC sampler from euroformix.
 - textEditor: GUI for text editor
 - profileSwapperGUI: GUI for profile selector
 - getModelSettings: Helpfunction for obtaining defined model settings.
 - setMarkerSettings: Helpfunction for obtaining per-marker settings.
 - getEnvirKit: Helpfunction for obtaining specified kit.
 - clusterOrder: Helpfunction for obtaining order of ranking samples by genotype similarity
 
 Fixed issues:
 - In MatchMatrix panel, values are now truncated when modifying table and "Truncate" option is selected.
 
#Changes version 1.7.3 (21.04.2020)
 Fixed bugs:
 - Create report crashed when trying to print reference profiles when no references were given.
 - It is now possible to import only reference profiles in a case.
 
#Changes version 1.7.2 (21.04.2020)
 - Both Mixtures and evidence empty MatchStatus is now always 'Compared'.
 - Renamed option in "Advanced options": 'References in MatchStatus can be shown in MatchList' (before this was about whether to compare single source profiles or not).
 - When references are removed from GUI they are also removed from MatchStatus and MatchList
 - The MatchNetwork will always show the match candidates in MatchList (also those in MatchStatus). No more dynamic update if thresholds are changed.
 - The Match Matrix has several sort options: Sort wrt column/row-names, or sort columns/rows wrt largest summed proportions over all samples.
 - It is now possible to "delete" evid samples under "Selected profile(s)". This can only be done before running 'Compare' for first time.
 - Expected PH plots are removed from report.

 - The user can now select a subset of samples by checking "Advanced -> Advanced options -> Profile selector"
 -> Import: The user can select which samples to import from the case (sometimes the user don't want to import all data in the case).
 -> Report: The user can select which samples to include in the report (sometimes the user don't want to report the results regarding some of the samples). 

 - A language module has been added: The R-function has been getLanguage added to the library. A excel spreadsheet is used to recognize words.
   -> In GUI Select language under Advanced -> Advanced options. 

 Fixed bugs:
 - When showing MatchNetwork before comparisons are carried out.
 - Avoids the crash when doing deconvolution both with and without stutter when AMEL is included in allele frequencies.
   -> Note: Using euroformix v3.0.0 makes it possible to do deconvolution for AMEL marker when stutter model is active.


#Changes version 1.6.1 (17.01.2020)
 - New feature "Calculate evidence concordance": Mix-mix allele comparison based on alleles. Proportion of coordance between all pairwise evid profiles. The MAC threshold is used to extract candidates. Results can be shown in report ("Show concordant evidence" added to report layout)
 - Single sources with corresponding assigned unknown will not be shown in "Matches" anymore (because it's already extracted from the corresponding profile).
 - Matches based only on allele comparison is now possible. LR can now be avoided by not defining pop frequencies (not recommended).
 - Minor changes:
  -> Previous calculated LR-results are removed when clicking compare (avoiding errors).
  -> Added row-index for IBS and evid. concordance comparisons.
  
#Changes version 1.5.0 (03.09.2020)
 - Backward compatibility: Older CaseSolver projects can now be opened in new version.
 - Single source profiles (not matching any refs) can now also be compared in "Compare" (by default). This can be turned on under "Advanced->Advanced options". 
 - CaseSolver now has a "SNP mode" which can be activated under "Advanced->Advanced options". 
  -> This will treat all evid profiles as mixtures, and number of contr is assumed equal 3 when using the Quan model.
 - CaseSolver now utilizes interactive plots using functions plotEPG2/plotMPS2/plotTopEPG2/plotTopMPS2 from the euroformix R-package (available from v2.2.0).
  -> Visualizing MPS data is now possible.
 - The MatchNetwork now also include an interactive plot (if plotly is installed)
 - In report: Evid profiles estimated as 1 contributors is now given under single source profiles.
 - The MatchMatrix is still calculated if popFreq not given (using loci from evid data)
 - GUI changes:
  -> In Data panel:
    - The "Functionalities" panel was moved above the tables. 
    - The button "Import Ref(s)" has been added: The user can import additional reference profiles to the case (requires a text file with genemapper format).
    - The button "Selected profiles" has been added: 
	-- Substituting the Deconvolution and Export button
	-- The selected/marked profiles in the tables can be selected for further operations: View/Export/Deconvolve/Delete
  -> In Deconvolution panel: The user can now export or remove selected candidates.
  -> In Matches panel: Dedicated "Show match networks" buttons were added to separate Mixture and Single source matches.
  -> "Mixtures" table now changed to "Matches".
  -> In Frequency file selection. The user can now deselect the frequency file ("Remove selected").

 - Minor changes:
   -- When closing window, the user is asked whether to quit instead of asked to save project.
   -- Create Report made more robust.
   -- A warning is given by user before running "Random IBS".
   -- New created project files becomes much smaller in size.
   -- Name changes: "Number of required optimizers (EFM)" instead of "Number of randoms in optimizer".

 - Bug fixes (thanks to Lourdes Prieto for finding these):
  -> In R-function tabToListRef: Reference allele names as NA caused problems.
  -> In R-function getMatchesLR in gui.R[Line 1528]: The matchlist argument of function calcQuanLRcomparison has to be a matrix.
 
#Changes version 1.4.1 (14.05.2019)
- MatchNetwork not updated after running quanLR in second round when modeltype=(qual,both) (thanks to Vibeke Bertelsen for finding the scenario).
-> An additional column named "type" was added to the resCompLR object (indicates whether the LR is based on qual/quan model).
- Fixed bug when performing deconvolution ("Deconvolve All") when having fewer number of contributors than number of conditional references (thanks to Lourdes Prieto for finding the scenario).  
- Small fixes in report when no match results are given.

#Changes version 1.4.0 (09.03.2019)
- Qual based LR made more robust: Optimizer now uses a start value out of 0.1,0.35,0.7 which maximizes likelihood function (instead of only 0.1).
- The user now has a choice whether to replace calculation when making a specific quan. LR calculations (up to user whether replacing LR or not). 
- The user can now choose whether the estimation of number of contributors should start with assuming one contributor (in "Advanced->Advanced options"). This is now default. Otherwise this will be the assigned number of contributors from the import.
- It is now possible to use project files between different users.
- Added sorting functionalities for different tables:
 -- The qualitative and quantitative Match lists can now be sorted wrt Evid and Ref names (in addition to LR values) 
 -- The Mixtures matchLists can now be sorted wrt Evid and Ref names (in addition to #id from data import) 
 -- Evidence profiles in the Data import table can now also be sorted wrt MatchStatus.
 -- The qualitative and quantitative matchLists can now be sorted wrt Evid and Ref names (in addition to LR values) 
- Minor bugs:
 -- The error caused when removing references such that only one remains is fixed. 
 -- The order of LRs for conditional contributors when evaluating with quan. model from Match List 1, was not always correct.

#Changes version 1.3.7 (19.12.2018)
  - DC conditioned on partial REFS now returns candidate shows deconvolved REF suggests
  - LR for conditional references from Qual.MatchList now also shown. The LRs are based on the method last used.

#Changes version 1.3.6 (22.11.2018)
  - BUG fixed in DC when conditoning on partial REFS.

#Changes version 1.3.5 (20.11.2018)
  - BUGs fixed when no refs are imported in case (report and add_ref)

#Changes version 1.3.4 (15.11.2018)
  - Homozygous genotypes for Y-markers are now given as only 1 allele for references.

#Changes version 1.3.3 (15.11.2018)
  - Bug in ncol(matchMATGUI[]) when having 1 ref (and 1 mix).
  - Bug in getStructuredData when having 2 alleles for Y-markers.

#Changes version 1.3.2 (14.11.2018)
 - Fixed bug in getStructuredData when having partial references by import.

#Changes version 1.3.1 (13.11.2018)
 - Faster structuring of data. Still limit in number of evidence (need to loop through all evidence).
   -- The software now handles MANY references (more than 1e5 depending on the memory).
     -> Handled in "edit/add references" in that user must select a segment of references.

 - The IBS calculations:
	-- Using only references matched in MatchStatus if more than 10000.
	-- Utilizes the R-package big.memory for extended memory handling (optional).
	-- Using only markers also found in popFreq if specified.

 - Bug(s) fixed:
   -- After running qualLR+quanLR, crashes when reanalysing comparison with negative qualLR not in quan matchlist.


#Changes version 1.2.3 (07.11.2018)

  #Added in report:
  - Option "Show Expected PH plots" has been added. If comparison is done, this will create EPG plots of the expected PH for all Mixtures with at least one matching reference (DC results also given if ran).

  #Added calculation/interpretation flexibilities:
  - Added button "Calculate Quan LRs" will caclulate the calculate the quantitative LR for all comparisons in the list. Results will come in Match list (Quan LR) - Final Mixture table will be updated wrt results.
  - From the list "Match list (Qual LR)", the user can now double click on a row to calculate Quan LR for a specific hypothesis (user specifies conditional references + numContr).
    -- Results will update the Match list (Quan LR) - Final Mixture table will be updated wrt results.
    -- The conditional references are now listed in the column after "numContr".
  - Added button in Mixtures "Deconvolve all samples" will automatically run DC for all samples.
  - "Deconvolution" added to Data panel. The user can now consider DC at any times.
   --  THe user can here select any number of replicates and condition of any refs as wanted.


  Minor updates:
  - The user is asked to restart the program when importing a new case if there were results from previous case.
  - The program asks user of saving project before crossing out window (Quit).
  - Graphical flexibility
    -- Added buttons in MatchMatrix (Change view, truncate values)
  - Y markers and other non-freq markers are no longer compared in match matrix (this earlier caused females to have small match comparisons).
   -- Now using loci specified in allele frequencies.
 - Under Setup->Model: Selecting model types now comes first.
 - Single alleles for Y-markers of unknowns from SS_evidence are not considered as homozygous anymore (checking if loc name starts with DY). 
   -- CODE: #if(av[1]==av[2] && (toupper(substr(loc,0,2))=="DY" || toupper(substr(loc,0,1))=="Y") ) av = av[1] 
 - "Add reference profile": The user can now edit the alleles of all references.
 - Added option in Report Layout:
   -- Optional whether to show "MatchStatus" column.   
 - Minor bug in Report fixed:  Mixture w/PH not shown when no SS given in data.


#Changes version 1.1.1 (18.09.2018)
- A bug in the function calcHp has been fixed.

#Changes version 1.1.0 (18.09.2018)
- The button "Change view" has been added to make it possible to view the allele data of evidence/reference vertically instead of horizontally 
- The user can now change model specification option for euroformix (degradation ON/OFF stutter ON/OFF)
  -> NB: AMEL locus is excluded if stutter=ON 
- The LR values for each comparison is labeled in the window "Deconvolution/Show expected peak heights". THe LR will be based on the last calculated one (Quan or Qual).

- The user can now use a frequency file with more loci systems than defined. Fixed by changing line 31 in calcQuanLRcomparison().

Thanks to Vibeke Bertelsen for discovering the following bugs:
- Bug when having only 1 match candidate in LR based comparison, when generating a report. Fixed by changing line 1299 and 1308 in gui().
- Bug when performing DC with locus drop-out: Fixed by changing line 843 in gui().


#Changes version 1.0.1 (03.08.2018)
- The bug "setupAdvanced could not be found", after saving "Advanced Setting", has been fixed.


#Historical versions:
 #v1.2: Added QuanLR and DC functionalities
 #v1.1: Small changes (gui,model)
 #v1.0: Community Release
 #v0.8: Matches either based on QualLR, QuanLR  or both 
 #v0.7: Supporting MPS kit (requires >=euroformix_1.11.0)
 #v0.5: Deconvolution added
 #v0.4: Calculate RMP added	
 #v0.3: Reports are added
 #v0.2: Tables are fixed
 #v0.1: Initial version

