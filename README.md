# MICB475_25W2_Team_10

Team 10: Mirren Buchanan, Gopika Makhija, Phoebe McNair-Luxon, Mia North, Madilyn Portas, Michaela Schweitzer
Zotero Group: https://www.zotero.org/groups/6428513/micb_475_t10 

## Team Meeting #1 February 2nd, 2026: Agenda

1. Discuss ideas for Project 2
   a. Phoebe's idea: Does host phylogeny constrain microbiome composition across fish species?
      What is the more important determinant of microbiome composition and diversity, the environment or host phylogeny?
       - This would use the fish dataset (from MICB 475)
   b. Michaela's idea: Does the predominant tree cover inhabiting BC affect the diversity of the microbial community when regenerating after OM removal?
       - This would use the soil dataset (from MICB 475)
   c. Anyone else's ideas 
2. Key to discuss: given our large group, are these ideas challenging enough?
3. How to flesh out the topic over the week/proposal preparation?
   a. Dividing up tasks
   b. Further research is needed to be happy with the topic 

## Team Meeting #1 February 2nd, 2026: Minutes

### Ideas for the project: 
- Fish dataset: There's one phylogenetic group that is overrepresented, everything else is sample size of 2. 
     - Looking at the ratios of morphology (how does the host trait affect microbiome composition?)
          - Would need to classify the fish according to these measures before analyzing the data
          - Would have to pick one body region (gill, gut, etc)
          - Would have to make the variables categorical (big fish, small fish, medium fish)
               - Could make the sample size really small
               - Every time you divide up the samples, you have to do more comparisons 
          - Separate into 4 datasets for each location of collection (then the question would become across which body site does morphology affect microbiome composition most
          - Maybe do a screening process of all the columns you think are interesting and do an alpha/beta diversity loop through all of them, see which gives you a significant p-value - then narrow the variables down to that (in 4 datasets: all the different body locations)
          - After the screening, do the diversity and functional analysis or building a machine learning model
          - There's a paper that did this process with the dysautonomia dataset (https://ojs.library.ubc.ca/index.php/UJEMI/article/view/199531/193234)
          - Could do Spearman's test to do correlations (keeping continuous variables)
          - Just run the modeling with alpha diversity, then run beta diversity when we come up with the variables of interest (at this point we would need categorical variables) 
- Soil dataset: A group already did tree cover unfortunately 

### Proposal prep notes:
- Need to process the dataset in QIIME2 BEFORE the proposal is due
- Then our work moves into R for the actual analysis
- Finding morphological dimension category to focus on for this project
- Proposal will be pretty vague because we don't know what's going to come out of the modeling (variable x, variable y)
- Why are we randomly studying fish without a specific goal in mind at this point? 

### Idea for project workflow:
- Phyloseq analysis in R (gives you Shannon diversity for everything)
- Aim 1: Divide into 4 datasets (gill, hind, mid, skin) - at this point it leaves all the species 
- Do a for loop to give you p-values for alpha diversity metrics (Spearman's rank or Kruskal-Wallis)
- This will hopefully identify something that's interesting
- Turn whatever variable into something categorical (big, medium, small)
- Aim 2: Run diversity metrics (beta!)
- Aim 3: Run indicator taxa or core microbiome or deseq
- Aim 4: Run functional analysis (picrust2) or random forest (depending on the category we find)

### For next week:
- Start getting our thoughts together for this
- Maybe starting coding
- Start framing the proposal
- Why are we doing each aim in relation to our whole exploratory research question?

## Team Meeting #2 February 9th, 2026: Agenda
1. Discuss preliminary data processing on QIIME2.
2. Discuss outline for key aims and rationales (to be expanded on in the proposal).
  - Link to proposal: https://ubcca-my.sharepoint.com/:w:/g/personal/micswtz_student_ubc_ca/IQDAbC8BS-N7Q7AE4KS3XcOfAc0gmxKv7MzBp9uRUK7lE9E?e=vWt2lS 
4. Discuss key highlights for the introduction (setting up the context for our project).
5. Outline and refine research question and hypothesis.
6. Quick discussion of lit review on previous UJEMI publication that followed a similar workflow as our planned project.

## Team Meeting #2 February 9th, 2026: Agenda
1. Discuss preliminary data processing on QIIME2.
2. Discuss outline for key aims and rationales (to be expanded on in the proposal).
3. Discuss key highlights for the introduction (setting up the context for our project).
4. Outline and refine research question and hypothesis.
5. Quick discussion of lit review on previous UJEMI publication that followed a similar workflow as our planned project.

## Team Meeting #2 February 9th, 2026: Minutes
### Data processing
- Showing Hans the demultiplexed graph:
  - Can't set our threshold at 30
  - Can't trim at 70 (would be way too low)
  - Setting at 218-220 for trimming
 
### Research question and Aims
- Regarding hypothesis - we shouldn't have a true hypothesis
  - Just say that you are exploring the variables
- Make sure in Aim 4 we're giving the option to do a functional analysis
- Research question: something to think about is using the word "predict"
  - We would need to have experiments to be able to test this "prediction" - it would HAVE to be machine learning (this way don't talk about functional)
  - If we wrote it as functional instead (for Aim 4) don't use the word predict (taxonomic AND functional diversity) 
  - Predict might not be best: we're looking for correlations (which variables impact diversity) - find a correlation between a specific variable and a specific change in diversity (taxonomic and functional) 
  - Take out the specific names of the body sites and just change generally to "body sites"
  - Change it to a "to what extent" question instead of "which"

