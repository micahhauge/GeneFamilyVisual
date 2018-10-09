# Meeting Notes

## July 24, 2018 - Updated To-dos

Development:
1. Backend code (pre-JSON)
	+ Complete development 
	+ Test
	+ Refine
2. Visualization (post-JSON)
	+ Update/Dynamic Title (Get from the user)
	+ Enlarge ovals on hover over (just like the circles)?
	+ About link
	+ Test
3. Create test data sets - Dr. Carroll
 	+ JSON format
		+ albuminGeneFamily
		+ albuminWithRandom
		+ randomGenes

Deployment:
1. Domain name registration
	+ Decide on a name - All
	+ Register - Dr. Carroll
2. Web server set-up - Micah
3. Web documentation (user's manual) - ?

Writing:
1. Write draft
	+ Outline
	+ Methods
		+ Backend
		+ Visualization
	+ Results
	+ Conclusion - Drs. Carroll & Seipelt-Thiemann [waiting on Methods & Results milestones]
	+ Introduction - Dr. Seipelt-Thiemann [waiting on Conclusion milestones]
	+ Abstract
2. Share with co-authors - Dr. Carroll
3. Revise draft
4. Return to 1.

Manuscript submission - Dr. Carroll


## Feb 15, 2018 - Updated To-dos

Development:
1. Backend code (pre-JSON) - Ethan [ ]
	+ Complete development 
	+ Test
	+ Refine
2. Visualization (post-JSON) - Micah [ ]
	+ Complete development
		+ Save as image option
			+ JPEG
		+ Automatically choose the colors if not present in JSON
	+ Test
	+ Refine (return to Test)
3. Create test data sets - Dr. Carroll
 	+ JSON format
		+ albuminGeneFamily
		+ albuminWithRandom
		+ randomGenes

Deployment:
1. Domain name registration
	+ Decide on a name - All
	+ Register - Dr. Carroll
2. Web server set-up - Micah
3. Web documentation (user's manual) - ?

Writing:
1. Write draft
	+ Methods
		+ Backend (Ethan) [ ]
		+ Visualization (Micah) [ ]
	+ Results (Ethan & Micah) [ ]
	+ Conclusion - Drs. Carroll & Seipelt-Thiemann [waiting on Methods & Results milestones]
	+ Introduction - Dr. Seipelt-Thiemann [waiting on Conclusion milestones]
	+ Abstract
2. Share with co-authors - Dr. Carroll
3. Revise draft
4. Return to 1.

Manuscript submission - Dr. Carroll


## Feb 5, 2018 - Updated To-dos

Development:
1. Backend code (pre-JSON) - Ethan [ ]
	+ Complete development 
	+ Test
	+ Refine
2. Visualization (post-JSON) - Micah [ ]
	+ Complete development
		+ Save as image option
			+ JPEG (and PNG)
			+ SVG
	+ Test
	+ Refine (return to Test)
3. Create test data sets - Dr. Carroll
 	+ Unrelated genes (just needs JSON format answer key)
	+ Add unrelated genes to a gene family (just needs JSON format answer key)

Deployment:
1. Domain name registration
	+ Decide on a name - All
	+ Register - Dr. Carroll
2. Web server set-up - Micah
3. Web documentation (user's manual) - ?

Writing:
1. Write draft
	+ Methods
		+ Backend (Ethan) [ ]
		+ Visualization (Micah) [ ]
	+ Results (Ethan & Micah) [ ]
	+ Conclusion - Drs. Carroll & Seipelt-Thiemann [waiting on Methods & Results milestones]
	+ Introduction - Dr. Seipelt-Thiemann [waiting on Conclusion milestones]
	+ Abstract
2. Share with co-authors - Dr. Carroll
3. Revise draft
4. Return to 1.

Manuscript submission - Dr. Carroll



## Nov 20, 2017 - MTSU Meetings Summary

### New JSON Record Format
Micah helped me realize today that the role of the visualization front-end is to just visualize and therefore to do as little computation as necessary.  Consequently, we need to change the JSON record format specification (from our March 20th notes) to be have the information already processed for layout.  Micah has agreed to lead out on that.  A rough sketch of what we talked about is:

List of visual exons
+ Label (length(s) for now)
+ Visual column index
+ Length (number of visual columns [most will be 1]
+ UTR type [0-2] (0: Not a UTR, 1: partial UTR, 2: full UTR)

List of gene info
+ For each gene:
  + Name
  + Visual exon indices

Micah's going to change his hard-coded example to use a JSON record.  (**Micah, can you have that done by the end of November 30th?**).

### Web server update
The prospects of using a server at CSU have not improved.  Leewood researched the pricing for web hosting companies.  We're looking for a Linux machine.  Additionally, we need to choose a domain name and register it.  (**Leewood, could you research the pricing and how to register a domain name before the end of the month?**)  (**All, who's got a great URL idea (that's not already taken)?**)


### Manuscript Outline
Dr. Seipelt-Thiemann and I had a great meeting this morning. We talked about the manuscript outline and made assignments for each section and sub-section.


### Milestones
Our goal is to submit our manuscript by the end of December.  To accomplish this, we'll need to meet a few more times.  I'll send out a Doodle poll to get everyone's availability.  Below are the remaining to-dos.  **If you see empty square brackets by your name, please fill-in a milestone goal for that item.**  As you complete an item, email everyone so that we can celebrate with you!


Development:
1. JSON record format - Micah [November 30, 2017]
2. Create test data sets - Dr. Carroll [December 7, 2017]
	+ Unrelated genes (Phase I done November 29, 2017)
	+ Add unrelated genes to a gene family (Phase I done November 29, 2017)
3. Backend code (pre-JSON) - Ethan [ ]
	+ Complete development 
	+ Test
	+ Refine (return to Test)
4. Visualization (post-JSON) - Micah [ ]
	+ Complete development
		+ Save as image option
			+ JPEG (and PNG)
			+ SVG
	+ Test
	+ Refine (return to Test)

Deployment:
1. Domain name registration
	+ Research pricing & how to - Leewood? [ ]
	+ Decide on a name - All
	+ Register - Dr. Carroll
2. Web server set-up - Dr. Carroll
3. Web documentation (user's manual) - ?

Writing:
1. ~~Venue selection - Dr. Carroll [December 1, 2017]~~ <Completed December 1, 2107>
2. Write draft
	+ Methods
		+ Backend (Ethan) [ ]
		+ Visualization (Micah) [ ]
	+ Results (Ethan & Micah) [ ]
	+ Conclusion - Drs. Carroll & Seipelt-Thiemann [waiting on Methods & Results milestones]
	+ Introduction - Dr. Seipelt-Thiemann [waiting on Conclusion milestones]
	+ Abstract
3. Share with co-authors - Dr. Carroll
4. Revise draft
5. Return to 2.

Manuscript submission - Dr. Carroll [December 30, 2017]


## Sep 14, 2017 - Google Hangout Notes

### Development:
1. Backend code (pre-JSON) - **Ethan**  
	A. Complete development  
	B. Create test data sets  
	C. Test  
	D. Refine (return to 1.C.)  
2. Visualization (post-JSON)  
	A. Complete development  
		i. Save as image option  
	B. Create test cases  
	C. Test  
	D. Refine (return to 2.C.)  

### Deployment:
1. Web server set-up  
	A. CSU, or  
	B. Independent site - **Leewood**  
2. Web documentation - **Leewood**  

### Writing:   
1. Outline paper - **Dr. Carroll**  
	A. Coordinate with Dr. Seipelt-Thiemann - **Dr. Carroll**  
2. Share outline with co-authors  
3. Revise outline  
4. Write draft (return to 2.)  

### Manuscript submission - **Dr. Carroll**  


### To-dos before next meeting:  
**Dr. Carroll**  
+ Outline paper (with Dr. Seipelt-Thiemann)  
+ Send out Doodle poll for bi-weekly meeting  
+ Keep waiting for response about CSU server to handle python CGI script  

**Ethan**:  
+ Test cases  


## Apr 17, 2017  

Next meeting at 3:00 PM  

### Status reports  
**Ethan**: Emailed: Close to fully exporting valid JSON records  
**Leewood**: Worked on a JSON demo  
**Micah**:   


### To-dos  
**Ethan**: Visual exon columns  
**Leewood**: Have a python CGI script make a new HTML page and redirect to it  
**Micah**: Make a ordered list of feature improvements  


## Apr 17, 2017  

### Status reports  
**Ethan**: Emailed  
**Leewood**: Emailed: Looked into sockets and files for transferring JSON string  
**Micah**: Worked on visualization  


### To-dos  
**Ethan**: Code to produce JSON records by April 24  
**Leewood**: Look into websockets and make a simple example  
**Micah**: Review what is uploaded to github; email **Dr. Carroll** (so that he can email Dr. Seipelt-Thiemann);   


## Apr 10, 2017  

### Status reports  
**Ethan**:   
**Leewood**: Emailed  
**Micah**:   


### To-dos  
**Ethan**: Code to produce JSON records by April 24  
**Leewood**:   
**Micah**:   



## Apr 3, 2017  

### Status reports  
**Leewood**: Emailed  
**Ethan**: Orange juice  
**Micah**: Algorithm for color selection  


### To-dos  
**Ethan**:   
**Micah**: Find module for color selection (e.g., given 20 colors, return 20 complimentary hex values), example(s) for visual exons spanning multiple columns (default to placing them at the top or bottom?)  
**Leewood**: Researching JSON format  

## Mar 27, 2017  

### Status reports  

### To-dos  
**Ethan**: Scholars' Week presentations; code to delimit visual exons  
**Micah**: Fine-tuning for presentations; Counter example(s) for visual exons spanning multiple columns  
**Leewood**: Researching JSON format  

## Mar 20, 2017  

JSON record format:  
---------------------------------  
A dictionary, with  
key0: "geneFamilyInfo", which is a dictionary, with gene names as keys, and a list of exon properties (dictionaries) as the values  
The exon properties dictionary with the following keys:  
+ start MSA alignment position  
+ end MSA alignment position  
+ visual exon index  
key1: "visualExons", which is a list with each element having the following data:  
	 + (ending) reading frame [0-2]  
	 + Visual exon columns (list of 2 elements, start and end [so most visual exons have the same start and end values])  
	 + UTR type [0-2] (0: Not a UTR, 1: partial UTR, 2: full UTR)  



## Mar 13, 2017  
JSON record format:  
---------------------------------  
dictionary, with gene names as keys, and a list of exon properties (dictionaries) as the values  
The exon properties dictionary with the following keys:  
+ Length  (number of nucleotides)  
+ (ending) reading frame [0-2]  
+ Visual exon columns (list of 2 elements, start and end [so most visual exons have the same start and end values])  
+ UTR type [0-2] (0: Not a UTR, 1: partial UTR, 2: full UTR)  
+ start MSA alignment position  
+ end MSA alignment position  
An additional key "exonSet", with a list of all unique exon properties (just reading frame and visual exon columns)  

The set of unique visual exons is governed by the reading frame and visual exon columns list  

### Status reports (emailed before 8 AM the morning of our meeting):  



## Feb 27, 2017  
### Status reports (emailed before 8 AM the morning of our meeting):  
Dr. C:     
**Arol**:    Emailed  
**Micah**:   Title, find & review JSON specs (from previous meeting notes)  
**Ethan**:   Flowchart (phase I)  
**Leewood**:   

### To-Dos:  
**Arol**:      
**Micah**:	 Test Cases for visual exons spanning multiple columns  
**Ethan**:   Update JSON example;  Refine the flowchart; contact Dr. C; work on implementing the flowchart  
**Leewood**:   



## Feb 13, 2017  
### Status reports (emailed before 8 AM the morning of our meeting):  
Dr. C:     
**Arol**:    CGI scripts  
**Micah**:   (submitted via email)  
**Ethan**:   Milestones (done), started flowchart  
**Leewood**:   

### To-Dos:  
**Arol**:      
**Micah**:   Title, visual exons spanning multiple columns, find & review JSON specs (from previous meeting notes)  
**Ethan**:   Flowchart  
**Leewood**: Comment CGI script  


## Feb 6, 2017  
### Status reports (emailed before 8 AM the morning of our meeting):  
Dr. C:     
**Arol**:    Background  
         Look at github  
**Micah**:   Milestones (in github)  
**Ethan**:   Milestones  
         + Flowchart of algorithm to determine which exons map to which visual exons  
**Leewood**:  

## 01/27/2017  

New Meeting time (that doesn't conflict with clubs?)  
+ Mondays, 2:30 PM  

### Status reports (due 8 AM the morning of our meeting):  
1. Tasks accomplished  
2. Stumbling blocks encountered  
3. Plans for the next week  

### To-Dos:  
Dr. C:     
**Arol**:      
+ Background  
+ Look at github  
**Micah**:   Milestones  
         Add all relevant files from **Leewood**'s github repo  
**Ethan**:   Milestones  
         + Flowchart of algorithm to determine which exons map to which visual exons  
