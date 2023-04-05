#include "Tools.h"

#if DEBUG
	#ifndef CK
		#define CK 1
	#endif

	#ifndef CK_2
		#define CK_2 0
	#endif
	
	#ifndef CK_3
		#define CK_3 0
	#endif
#endif

#ifndef PRINT_TABLE
	#define PRINT_TABLE 0
#endif

//Input parameters
dataTypelenSeq k_min=16;
double tau=1.0/2.0;

dataTypeNChar n=0; //eBWT size

#define VALUE_F 100

              //A B C D E F G H I J K L M N O P Q R S T U V W X Y
int ORD[25] = {0,0,1,0,0,0,2,0,0,0,0,0,0,4,0,0,0,0,0,3,0,0,0,0,0};
#define ord(c) (ORD[c-65])
//terminator character at the end of the reads
char TERM = '#';

dataTypeNSeq numData=0;

vector<string> NAME;

bool Close(dataTypeNChar ind_start, dataTypeNChar ind_fin, FILE *inCDA,FILE *inEBWT,type_node &isIn, dataTypeNChar *n_clust){
	bool res = false;
	//ind_start: starting index of the eBWT positional cluster
	//ind_fin: final index (NOT included)
	//isIn stores colors in the cluster
	//B_Symb stores symbols in the cluster
	std::bitset<ALF> B_Symb;
	
	//To upload in memory positional cluster in eBWT and cda
	dataTypeNChar c_len=ind_fin-ind_start;
	dataTypeNSeq *bufferCDA = new dataTypeNSeq[c_len];
	dataTypedimAlpha *buffereBWT = new dataTypedimAlpha[c_len];
	
	fseek(inCDA,ind_start*sizeof(dataTypeNSeq), SEEK_SET);
	dataTypeNChar numchar=fread(bufferCDA,sizeof(dataTypeNSeq),c_len,inCDA);
	assert(numchar==c_len);
	fseek(inEBWT,ind_start*sizeof(dataTypedimAlpha), SEEK_SET);
	numchar=fread(buffereBWT,sizeof(dataTypedimAlpha),c_len,inEBWT);
	assert(numchar==c_len);
	
	for(dataTypeNChar index=0; index<c_len; index++){
		isIn[bufferCDA[index]]=1;
		if(buffereBWT[index]==TERM){
			B_Symb[ALF-1]=1;
		}
		else{
			B_Symb[ord(buffereBWT[index])]=1;
		}
	}
	
	//Check positional cluster maximal
	if(B_Symb.count()>1){
		#if CK_3 
			fprintf(stderr,"Close returns 1 -->\nCDA:\t");				
			for(dataTypeNChar index=ind_start; index<ind_fin; index++)
				fprintf(stderr,"%d", (int)bufferCDA[index]);
			fprintf(stderr,"\neBWT:\t");
			for(dataTypeNChar index=ind_start; index<ind_fin; index++)
				fprintf(stderr,"%c",buffereBWT[index]);
			fprintf(stderr,"\n");
		#endif
		(*n_clust)++;
		res=true;
	}
	
	delete [] bufferCDA;
	delete [] buffereBWT;
	
	return res;
}

bool CloseAndReduce(dataTypeNChar ind_start, dataTypeNChar ind_fin, FILE *inCDA,FILE *inEBWT, vector<dataTypeNSeq> &colors,type_node &isIn, dataTypeNChar *n_clust){
	bool res = false;
	//ind_start: starting index of the eBWT positional cluster
	//ind_fin: final index (NOT included)
	//isIn stores colors in the cluster
	//B_Symb stores symbols in the cluster
	std::bitset<ALF> B_Symb;
	
	std::pair<std::vector<dataTypeNSeq>::iterator,std::vector<dataTypeNSeq>::iterator> color_bounds; 
		
	//To upload in memory positional cluster in eBWT and cda
	dataTypeNChar c_len=ind_fin-ind_start;
	dataTypeNSeq *bufferCDA = new dataTypeNSeq[c_len];
	dataTypedimAlpha *buffereBWT = new dataTypedimAlpha[c_len];
	
	fseek(inCDA,ind_start*sizeof(dataTypeNSeq), SEEK_SET);
	dataTypeNChar numchar=fread(bufferCDA,sizeof(dataTypeNSeq),c_len,inCDA);
	assert(numchar==c_len);
	fseek(inEBWT,ind_start*sizeof(dataTypedimAlpha), SEEK_SET);
	numchar=fread(buffereBWT,sizeof(dataTypedimAlpha),c_len,inEBWT);
	assert(numchar==c_len);
	
	for(dataTypeNChar index=0; index<c_len; index++){
		color_bounds=std::equal_range (colors.begin(), colors.end(), bufferCDA[index]); 
		if(color_bounds.second-color_bounds.first==1){
			isIn[bufferCDA[index]]=1;
			if(buffereBWT[index]==TERM){
				B_Symb[ALF-1]=1;
			}
			else{
				B_Symb[ord(buffereBWT[index])]=1;
			}
		}
	}
	
	//Check positional cluster maximal
	if(B_Symb.count()>1){
		#if CK_3 
			fprintf(stderr,"Close returns 1 -->\nCDA:\t");				
			for(dataTypeNChar index=ind_start; index<ind_fin; index++)
				fprintf(stderr,"%d", (int)bufferCDA[index]);
			fprintf(stderr,"\neBWT:\t");
			for(dataTypeNChar index=ind_start; index<ind_fin; index++)
				fprintf(stderr,"%c",buffereBWT[index]);
			fprintf(stderr,"\n");
		#endif
		(*n_clust)++;
		res=true;
	}
	
	delete [] bufferCDA;
	delete [] buffereBWT;
	
	return res;
}

int UpdateTable(std::unordered_map<bitset<SIZE_BITSET>,dataTypeNChar> &votes, dataTypelenSeq w, type_node &isIn, nodes &inSet){
	//w is the weight to use for updating partial_sum
	#if CK
		fprintf(stderr,"\nUpdateTable --> isIn: ");				
		for(dataTypeNSeq s=0; s<inSet.size(); s++)
			fprintf(stderr,"%d,", (int)isIn[s]);
		fprintf(stderr,"\n-->c_weight=%d\n",w);
	#endif
	
	//B_vect stores which elements of inSet are in the cluster
	type_node B_vect;
	
	//the cluster votes the candidate part if ALL the elements appear over the threshold tau
	bool vote=true;
	double abund=0.0;

	//scan inSet
	dataTypeNSeq k=0;
	while ( vote && k<inSet.size()){
		dataTypeNSeq card=inSet[k].count();
		dataTypeNSeq i=(isIn&inSet[k]).count();
		abund=(double)i/(double)card;
		#if CK_3
			fprintf(stderr,"\t--> inSet[%d]: ",k);				
			for(dataTypeNSeq s=0; s<numData; s++)
				fprintf(stderr,"%d,", (int)inSet[k][s]);
			fprintf(stderr,"\n-->abund=%lf\n",abund);
		#endif
		if(abund>=tau)
			B_vect[k]=1;
		else if(abund>0)
			vote=false;

		k++;
	}
	#if CK
		fprintf(stderr,"vote -->%d\n",vote);
	#endif
	
	//Increase table
	if(vote && (B_vect.count()<inSet.size()) && (B_vect.count()>1)) {
		dataTypeNChar index = B_vect.to_ulong();
		
		auto search_ele = votes.find(index);
		if(search_ele != votes.end()){
			#if CK
			fprintf(stderr,"Update votes for index = %lu --> search_ele->first: %lu, search_ele-> second: %lu\n",index,search_ele->first,search_ele->second);				
			#endif
			search_ele->second = search_ele->second + w;
			#if CK
			fprintf(stderr,"Update votes --> search_ele->first: %lu, search_ele-> second: %lu\n",search_ele->first,search_ele->second);				
			#endif
		}
		else{
			votes.insert({index,w});
			#if CK
			fprintf(stderr,"New entry votes for index = %lu --> votes[index]: %lu\n",index,votes[index]);				
			#endif
		}
		return 1;
	}
	else return 0;
	
}

//when this is over, votes_list contains a list of entries that each vote 1, each entry is a group of elements of ele_in_set
dataTypeNChar ClusterAnalysis(std::unordered_map<bitset<SIZE_BITSET>,dataTypeNChar> &votes,nodes &ele_in_set,FILE* fileInLCP,FILE *fileInCDA_2,FILE *fileInCDA,FILE *fileInEBWT, dataTypeNChar *n_clust){
	
	dataTypeNChar nClusters=0;
	
	//colors is a vector containing the IDs of CDA elements to analyze
	//if (colors.size()<numData) DS must be reduced
	bool reduce = false;
	
	//Ordered vector of CDA elements in ele_in_set
	vector<dataTypeNSeq> colors;
	type_node union_in_set;
	union_in_set.reset();
	for(dataTypeNSeq i=0; i<ele_in_set.size(); i++){
		union_in_set|=ele_in_set[i];
	}
	for(dataTypeNSeq i=0; i<numData; i++){
		if(union_in_set[i]==1) 
			colors.push_back(i);
	}
	
	if(colors.size()<numData)
		reduce=true;
	
	#if CK
		fprintf(stderr,"ClusterAnalysis --> ele_in_set.size=%lu, colors=[",ele_in_set.size());
		for(dataTypeNSeq i=0; i<colors.size(); i++)
			fprintf(stderr,"%d,", (int)colors[i]);
		fprintf(stderr,"]\n");
	#endif
	
	bool started=false;
	dataTypelenSeq pred = -1;
	bool grow = true;
	dataTypeNChar tail_len = 0;
	dataTypeNChar pos_init = 0, pos_end = 0;
	dataTypelenSeq c_weight=1;

	//Scan LCP file to detect maximal eBWT positional clusters
	//To read LCP and CDA
	dataTypeNChar numcharLCP,numcharCDA;
	dataTypeNChar index_file=0;
	dataTypelenSeq *bufferLCP = new dataTypelenSeq[BUFFERSIZE];
	dataTypeNSeq *bufferCDA = new dataTypeNSeq[BUFFERSIZE];
	
	fseek(fileInLCP, 0, SEEK_SET);
	
	//Read BUFFERSIZE elements from LCP-array
	numcharLCP=fread(bufferLCP,sizeof(dataTypelenSeq),BUFFERSIZE,fileInLCP);
	
	if(not reduce){
		while(numcharLCP>0){
			for(dataTypeNChar indexbuffer=0; indexbuffer<numcharLCP; indexbuffer++){	
				if(bufferLCP[indexbuffer]>=k_min) {
					if (not started){
						started=true, grow=true;
						pos_init=index_file+indexbuffer-1;
						tail_len=0;
					}
					else //cluster is open
					{
						//Case1: increasing sequence of LCP values
						if (grow){
							if(bufferLCP[indexbuffer]<pred) grow = false;
						}
						//Case2:non-increasing sequence of LCP values
						else{ 
							if (bufferLCP[indexbuffer]<pred) tail_len=0; 
							else if (bufferLCP[indexbuffer]==pred) tail_len++;
							else { //close cluster and start a new one
								type_node isIn;
								pos_end = index_file + indexbuffer-1-tail_len;
								if(Close(pos_init,pos_end,fileInCDA,fileInEBWT,isIn,n_clust)==1){	
									//Close returns 1 if the eBWT positional cluster is maximal					
									//Update votes
									//isIn: all colors in the cluster (not divided)
									nClusters+=UpdateTable(votes,c_weight,isIn,ele_in_set); 
								}
								//Start a new cluster
								started=true, grow=true;
								pos_init=pos_end;
								tail_len=0;
							}
						}
					}//end-else if(not started)		
				}
				else if (started){ //close cluster
					type_node isIn;
					pos_end = index_file + indexbuffer;
					if(Close(pos_init,pos_end,fileInCDA,fileInEBWT,isIn,n_clust)==1){
						nClusters+=UpdateTable(votes,c_weight,isIn,ele_in_set);
					}
					started=false;
				}//end-else if(b_lcp[indexbuffer]>=k_min)
				pred=bufferLCP[indexbuffer];
			}//end-for
		
			index_file+=numcharLCP;
			//Read BUFFERSIZE elements from LCP-array and CDA-array
			numcharLCP=fread(bufferLCP,sizeof(dataTypelenSeq),BUFFERSIZE,fileInLCP);
			numcharCDA=fread(bufferCDA,sizeof(dataTypeNSeq),BUFFERSIZE,fileInCDA_2);
		}//end-while
			
		//Close a possibly last cluster 
		if(started){
			type_node isIn;
			if(Close(pos_init,index_file,fileInCDA,fileInEBWT,isIn,n_clust)==1){
				nClusters+=UpdateTable(votes,c_weight,isIn,ele_in_set);
			}
		}	
	}//end-if not reduce
	else{
		//To restrict DS
		dataTypelenSeq inherited_lcp=MAX_LCP;
		bool prev_removed=false;
		std::pair<std::vector<dataTypeNSeq>::iterator,std::vector<dataTypeNSeq>::iterator> color_bounds;
	
		dataTypeNChar pred_index = 0;
		
		fseek(fileInCDA_2, 0, SEEK_SET);
		numcharCDA=fread(bufferCDA,sizeof(dataTypeNSeq),BUFFERSIZE,fileInCDA_2);
		assert(numcharLCP==numcharCDA);
	
		while(numcharLCP>0){
			for(dataTypeNChar indexbuffer=0; indexbuffer<numcharLCP; indexbuffer++){	
				color_bounds=std::equal_range (colors.begin(), colors.end(), bufferCDA[indexbuffer]);  
				//if bufferCDA[indexbuffer] IS present--> check LCP entry and analize it
				if(color_bounds.second-color_bounds.first==1){
					//Possibly change the LCP entry
					if(prev_removed)	
						bufferLCP[indexbuffer]=(bufferLCP[indexbuffer]<inherited_lcp)?bufferLCP[indexbuffer]:inherited_lcp;
					prev_removed=false;
					inherited_lcp=MAX_LCP;
					//Analize the LCP entry
					if(bufferLCP[indexbuffer]>=k_min) {
						if (not started){
							started=true, grow=true;
							pos_init=pred_index;
							tail_len=0;
						}
						else //cluster is open
						{
							//Case1: increasing sequence of LCP values
							if (grow){
								if(bufferLCP[indexbuffer]<pred) grow = false;
							}
							//Case2:non-increasing sequence of LCP values
							else{ 
								if (bufferLCP[indexbuffer]<pred) tail_len=0; 
								else if (bufferLCP[indexbuffer]==pred) tail_len++;
								else { //close cluster and start a new one
									type_node isIn;
									if(CloseAndReduce(pos_init,pred_index,fileInCDA,fileInEBWT,colors,isIn,n_clust)==1){	
										//Close returns 1 if the eBWT positional cluster is maximal					
										//Update votes
										//isIn: all colors in the cluster (not divided)
										nClusters+=UpdateTable(votes,c_weight,isIn,ele_in_set); 
									}
									//Start a new cluster
									started=true, grow=true;
									pos_init=pred_index;
									tail_len=0;
								}
							}
						}//end-else if(not started)		
					}
					else if (started){ //close cluster
						type_node isIn;
						pos_end = index_file + indexbuffer;
						if(CloseAndReduce(pos_init,pos_end,fileInCDA,fileInEBWT,colors,isIn,n_clust)==1){
							nClusters+=UpdateTable(votes,c_weight,isIn,ele_in_set);
						}
						started=false;
						tail_len=0;
					}//end-else if(bufferLCP[indexbuffer]>=k_min)
					
					pred=bufferLCP[indexbuffer];
					if(tail_len==0)
						pred_index=index_file+indexbuffer;
				}
				//if bufferCDA[indexbuffer] is NOT --> do NOT write, but SAVE inherited_lcp
				else{
					prev_removed=true;
					inherited_lcp=(bufferLCP[indexbuffer]<inherited_lcp)?bufferLCP[indexbuffer]:inherited_lcp;
				}
			}//end-for
			
			index_file+=numcharLCP;
			//Read BUFFERSIZE elements from LCP-array and CDA-array
			numcharLCP=fread(bufferLCP,sizeof(dataTypelenSeq),BUFFERSIZE,fileInLCP);
			numcharCDA=fread(bufferCDA,sizeof(dataTypeNSeq),BUFFERSIZE,fileInCDA_2);
		}//end-while
		
		//Close a possibly last cluster 
		if(started){
			type_node isIn;
			if(CloseAndReduce(pos_init,index_file,fileInCDA,fileInEBWT,colors,isIn,n_clust)==1){
				nClusters+=UpdateTable(votes,c_weight,isIn,ele_in_set);
			}
		}	
	}//end-else not reduce
	
	assert(index_file==n);
	
	delete [] bufferLCP;
	delete [] bufferCDA;
	
	return nClusters;
}

type_node inset_to_colors(type_node &curr, nodes &in_set){

	//Build the bitset set_ele corresponding to curr
	type_node set_ele;
	set_ele.reset();
	for(dataTypeNSeq s=0; s<in_set.size(); s++){
		if(curr[s]==1){
			//Add to set_ele in_set[s]
			set_ele|=in_set[s];
		}
	}
	return set_ele;
}

int insert_new_ele(tree<type_node>& tr, tree<type_node>::iterator iRoot, type_node new_ele)
{
	vector<tree<type_node>::iterator> subtrees_list;
	//Scan child1, ..., childn
	tree<type_node>::sibling_iterator iChildren;
	for (iChildren = tr.begin(iRoot); iChildren != tr.end(iRoot); ++iChildren) 
	{
		//Check intersection of new_ele and *iChildren
		type_node res = (new_ele&(*iChildren));
		
		if(res.any())//non-empty intersection
		{
			//Case 1: *iChildren contained in new_ele
			if(res.count()==(*iChildren).count()) 
				subtrees_list.push_back(iChildren);
			//Case 2: new_ele contained in *iChildren (N.B. there exists only one child with such a property)
			else if(res.count()==new_ele.count()) 
				return insert_new_ele(tr,iChildren,new_ele);
			//Case 3: *iChildren has ONLY some elements of new_ele --> we can discard the new_ele (it will not be a node of our tree)
			else
				return 0;
		}
	}
	
	//If only Case 1 has occurred
	int i, siblingNum=subtrees_list.size();
	
	//Check subtrees_list has more than 1 element (by construction)
	if(siblingNum<2)
		cerr << "EXCEPTION subtrees_list.size: " << subtrees_list.size() << endl;
	
	#if CK_2
	tr.debug_verify_consistency();
	#endif
	
	//TREE UPDATING...
	//Create a new tree having new_ele as root
	tree<type_node> new_tr;
	tree<type_node>::iterator it = new_tr.set_head(new_ele);
	//Add to new_tr children/Remove nodes from tr
	for(i=0; i<siblingNum; ++i){
		//Remove subtree rooted at subtrees_list[i] using move_out
		auto sub_tr = tr.move_out(subtrees_list[i]);
		
		#if CK_2
		sub_tr.debug_verify_consistency();	
		tr.debug_verify_consistency();
		#endif
		
		//Add as i-th child sub_tr to new_tr
		new_tr.move_in_as_nth_child(it, i, sub_tr);
		
		#if CK_2
		new_tr.debug_verify_consistency();
		cout << "new_tr:\n";
		kptree::print_tree_bracketed(new_tr);
		cout << "\ntr:\n";
		kptree::print_tree_bracketed(tr);
		cout << "\n";
		#endif
	}
	//Add new_tr as last child of iRoot
	tr.move_in_below(iRoot,new_tr);
	
	#if CK_2
	tr.debug_verify_consistency();
	//Print
	kptree::print_tree_bracketed(tr);
	cout << "\n";
	#endif
	
	return 1;
}
#if PRINT_TABLE
dataTypeNChar Update_phylotree(nodes &in_set,vector<type_entry> &list_subsets,tree<type_node>& phy_tr,tree<type_node>::iterator c_root, FILE *f_out)
#else
dataTypeNChar Update_phylotree(nodes &in_set,vector<type_entry> &list_subsets,tree<type_node>& phy_tr,tree<type_node>::iterator c_root)
#endif
{
	dataTypeNChar num_it = 0;
	dataTypeNSeq num_success = 0;
	dataTypeNChar consecutive_fail = 0;
	dataTypeNSeq n_max=in_set.size()-1;
	dataTypeNChar value_f = (VALUE_F>2*n_max)?2*n_max:VALUE_F;
	double max_value = log10(list_subsets[list_subsets.size()-1].second);
	
	bool stop=false;
	#if PRINT_TABLE
	fprintf(f_out,"PrintTable\n");
	#endif
	
	while((!stop) && (num_it<list_subsets.size())){
		//bit_max_ele is the maximum combination set among the remaining
		type_node bit_max_ele = inset_to_colors(list_subsets[list_subsets.size()-1-num_it].first,in_set);
			
		#if CK_2
			fprintf(stderr,"findMax --> num_it=%lu, bit_max_ele=",num_it);
			for(dataTypeNSeq s=0; s<numData; s++)
				fprintf(stderr,"%d,",(int)bit_max_ele[s]);
		#endif
		
		if(insert_new_ele(phy_tr,c_root,bit_max_ele)==1){
			#if CK_2
			fprintf(stderr,"Success - bit_max_ele inserted\n");
			#endif
			num_success++;
			consecutive_fail=0;
			#if PRINT_TABLE
			for(dataTypeNSeq s=0; s<numData; s++){
				if(bit_max_ele[s]==1){
					//Print the sth NAME
					fprintf(f_out,"%s,", &NAME[s][0]);
				}
			}
			fprintf(f_out,"\t%lu\t*\n",list_subsets[list_subsets.size()-1-num_it].second);
			#endif
		
			
		}
		else{
			#if CK_2
			fprintf(stderr,"Failure - bit_max_ele NOT inserted\n");
			#endif
			consecutive_fail++;
			#if PRINT_TABLE
			for(dataTypeNSeq s=0; s<numData; s++){
				if(bit_max_ele[s]==1){
					//Print the sth NAME
					fprintf(f_out,"%s,", &NAME[s][0]);
				}
			}
			fprintf(f_out,"\t%lu\t\n",list_subsets[list_subsets.size()-1-num_it].second);
			#endif
		}
		
		num_it++;
		
		if((consecutive_fail>=value_f) || (num_success>=n_max-1) || (max_value - log10(list_subsets[list_subsets.size()-1-num_it].second) > 2))
			stop=true;
	}//end-while
	
	fprintf(stderr,"END Update_phylotree --> N. iterations = %lu, N. Success = %u\n",num_it,num_success);
	
	return num_success;
}

#if PRINT_TABLE
bool Refinement(nodes &in_set,tree<type_node>::iterator root,tree<type_node>& tr,FILE* InFileLCP,FILE *InFileCDA, FILE *InFileCDA2,FILE *InFileEBWT, FILE *f_out)
#else
bool Refinement(nodes &in_set,tree<type_node>::iterator root,tree<type_node>& tr,FILE* InFileLCP, FILE *InFileCDA,FILE *InFileCDA2,FILE *InFileEBWT)
#endif
{
	time_t t_refine=0, t_total=0;
	clock_t c_refine=0, c_total=0;
	
	time_start(&t_refine, &c_refine); //start time
	time_start(&t_total, &c_total); //start time
	
	//map to store scores of the combination sets
	std::unordered_map<bitset<SIZE_BITSET>,dataTypeNChar> votes;
	
	//Detect positional clusters and write Table partial_sum
	dataTypeNChar n_clust=0;
 
	dataTypeNChar n_clust_vote=ClusterAnalysis(votes,in_set,InFileLCP,InFileCDA,InFileCDA2,InFileEBWT,&n_clust); 
	
	fprintf(stderr,"Refinement -> Number maximal positional clusters=%lu, of which voting=%lu\n",n_clust,n_clust_vote);
	//fprintf(stdout,"Refinement -> Number maximal positional clusters=%lu, of which voting=%lu\n",n_clust,n_clust_vote);
	//fprintf(stdout,"Refinement -> votes_list size: %lu\n", votes.size());
	//fprintf(stdout,"Time: %.6lf\n", time_stop(t_refine, c_refine));
	
	dataTypeNChar num_inserted = 0;
	
	if(votes.size()>0){
		//Order votes by score
		
		// Declare vector of pairs
		vector<type_entry> table;
		//pair<bitset<SIZE_BITSET>,dataTypeNChar> is type_entry
		//pair<set,score>
		
		// Copy key-value pair from Map
		for (auto& it : votes) {
			table.push_back(it);
		}
		//Clear votes
		votes.clear();
		
		// Sort using cmp_function
		sort(table.begin(), table.end(), cmp_function);

		//fprintf(stdout,"Refinement -> table size: %lu\n", table.size());
		//fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
		
		//Take the top combination sets to update the subtree of phylotree rooted at root	
		#if PRINT_TABLE
		num_inserted = Update_phylotree(in_set,table,tr,root,f_out);
		#else
		num_inserted = Update_phylotree(in_set,table,tr,root);
		#endif
		
		//Clear table
		table.clear();
	}
	
	#if CK_2
	fprintf(stderr,"END Refinement --> num_inserted=%lu\n",num_inserted);
	#endif
	
	return (num_inserted>0);
}
bool node_in_list(tree<type_node>::iterator pn, nodes &in_set){
	for(dataTypeNSeq i=0; i<in_set.size(); i++){
		//Compare *pn to in_set[i]
		if(*pn==in_set[i])
			return true;
	}
	return false;
}
void Update_job_list(std::deque<job>& job_list,tree<type_node>::iterator root,tree<type_node>& tr, nodes &in_set)
{
	#if CK_2
		cout << "Visit tree node: " << *root << endl;
		cout << "N. children: " << tr.number_of_children(root) << endl;
	#endif
		
	if (tr.number_of_children(root)==0 || node_in_list(root,in_set)) //Stop recursion
		return;
	else
	{
		//Recursively call Update_job_list
		//In case tr.number_of_children(root) > 2: append root and its children to job_list
		nodes upperSubTree;
		//Scan child1, ..., childn
		tree<type_node>::sibling_iterator child;
		for (child = tr.begin(root); child != tr.end(root); ++child) {
			//Recursive call
			Update_job_list(job_list,child,tr,in_set);
			if (tr.number_of_children(root) > 2)
				upperSubTree.push_back(*child);
		}
		if (tr.number_of_children(root) > 2){
			//Put the set of children in a new pair (upperSubTree,root) and append it to job_list
			job_list.push_back(make_pair(upperSubTree,root));
		}
	}
}

void Print_tree_newick_form(tree<type_node>& t,tree<type_node>::iterator root, std::ostream& str) 
{
	if(t.empty()) return;
	if (t.number_of_children(root) == 0){
		//Case leaf --> print its name
		dataTypeNSeq i=0;
		while(i<numData){
			if((*root)[i]){
				str << NAME[i];
				i=numData;
			}
			else
				i++;
		}
	}
	else {
		// parent
		str << "(";
		// child1, ..., childn
		dataTypeNSeq siblingCount = t.number_of_siblings(t.begin(root));
		dataTypeNSeq siblingNum;
		tree<type_node>::sibling_iterator children;
		for (children = t.begin(root), siblingNum = 0; children != t.end(root); ++children, ++siblingNum) {
			// recursively print child
			Print_tree_newick_form(t,children,str);
			// comma after every child except the last one
			if (siblingNum != siblingCount )
				str << ",";
		}
		str << ")";
	}
}

int main(int argc, char **argv) {
	//time_t t_refine=0, t_total=0;
	time_t t_total=0;
    //clock_t c_refine=0, c_total=0;
    clock_t c_total=0;
	
	//INPUT
	if( argc < 6) {
		fprintf(stderr,"Error usage %s fileFasta fileInfo fileOutput k_min tau\n",argv[0]); 
		exit(1);
	}

	string fileFasta=argv[1];
	string fileInfo=argv[2];
	string output_new=argv[3];

	sscanf(argv[4], "%hhu", &k_min);
	sscanf(argv[5], "%lf", &tau);
	
	//print input parameters
	fprintf(stdout,"k_min: %d\ntau: %lf\n",k_min,tau);
	
	fprintf(stderr,"sizeof(dataTypeNChar)=%lu bytes.\n",sizeof(dataTypeNChar));
	fprintf(stderr,"sizeof(dataTypelenSeq)=%lu bytes.\n",sizeof(dataTypelenSeq));
	fprintf(stderr,"sizeof(dataTypeNSeq)=%lu bytes.\n",sizeof(dataTypeNSeq));
	
	//Read fileInfo
	std::ifstream in_list;
	in_list.open(fileInfo.c_str(), std::ifstream::in);
	if (!in_list.is_open()){
		fprintf(stderr,"Error opening file %s\n",fileInfo.c_str());
		exit (EXIT_FAILURE);
	}
	
	while(getline(in_list,fileInfo,'\t')){
		NAME.push_back(fileInfo);
		getline(in_list,fileInfo,'\n');
		numData++;
	}
	
	fprintf(stdout,"numData: %u\n",numData);
	
	if(numData>SIZE_BITSET){
		fprintf(stderr,"ERROR numData is %u and SIZE_BITSET is %u\nPlease increase SIZE_BITSET in Tools.h",numData,SIZE_BITSET);
		exit(EXIT_FAILURE);
	}
	
    //Files .lcp, .ebwt and .cda 
	string fnLCP=fileFasta+".lcp";
	string fnCDA=fileFasta+".cda";
	string fnEBWT=fileFasta+".ebwt";
	
	FILE *InLCP= fopen(fnLCP.c_str(), "rb");
	FILE *InCDA= fopen(fnCDA.c_str(), "rb");
	FILE *InCDA_2= fopen(fnCDA.c_str(), "rb");
	FILE *InEBWT= fopen(fnEBWT.c_str(), "rb");
	
    if(InLCP==NULL){
		fprintf(stderr,"Error opening file %s\n",fnLCP.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    if(InCDA==NULL){
        fprintf(stderr,"Error opening file %s\n",fnCDA.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
	if(InCDA_2==NULL){
        fprintf(stderr,"Error opening file %s\n",fnCDA.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    if (InEBWT==NULL){
        fprintf(stderr,"Error opening file %s\n",fnEBWT.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
	}
	
	//Set the eBWT length
    fseek(InEBWT, 0, SEEK_END);
	n = ftell(InEBWT)/sizeof(dataTypedimAlpha);
	
	//Output file(s)
	filebuf out_fb;
	out_fb.open (output_new.c_str(),std::ios::out);
	ostream OutNew(&out_fb);
	//FILE *OutNew= fopen(output_new.c_str(), "w");
    #if PRINT_TABLE
	std::stringstream ssout;
	ssout << "Table_k_" << (int)k_min << "_tau" << (int)(tau*10) << ".txt\0";
	
	string fnOUT=ssout.str();
	FILE *OutP= fopen(fnOUT.c_str(), "w");
    if ((OutP==NULL)){
		std::cerr << "Error opening " << fnOUT << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
    #endif
    
	time_start(&t_total, &c_total); //start time

	//Queue of subsets of colors that must be partitioned
	std::deque<job> job_list;
	
	//Starting tree
	tree<type_node> phylo_tree;
	
	type_node ele; //root ele
	for(dataTypeNSeq i=0; i<numData; i++)
		ele[i]=1;

	//Initialize phylo_tree
	tree<type_node>::iterator root = phylo_tree.insert(phylo_tree.begin(), ele);
	
	//At the beginning the queue job_list has only one element --> non-final root
	nodes startSet;
	for(dataTypeNSeq i=0; i<numData; i++){
		type_node o;
		o[i]=1;
		startSet.push_back(o);
		//Initialize root children
		phylo_tree.append_child(root, o);
	}
	job_list.push_back(make_pair(startSet,root));
	
	#if CK_2
	//Check phylo_tree initialization
	kptree::print_tree_bracketed(phylo_tree, std::cout);
	std::cout << std::endl;
	#endif
	
	dataTypeNSeq step=0;
	
	while(job_list.size()>0){
		//Pick the first element and remove it from the queue
		job P = job_list.front();
		job_list.pop_front();
		
		//P.first is the set of nodes to merge together, P.second is a pointer to their parent in phylo_tree
		assert((P.first).size()>2);
		
		//Call Refinement
		//Refinement returns true if some nodes have been added to phylo_tree
		#if PRINT_TABLE
		if(Refinement(P.first,P.second,phylo_tree,InLCP,InCDA,InCDA_2,InEBWT,OutP))
		#else
        if(Refinement(P.first,P.second,phylo_tree,InLCP,InCDA,InCDA_2,InEBWT))
        #endif
			//Load new jobs
			Update_job_list(job_list,P.second,phylo_tree,P.first);
		//Otherwise nothing else is done -> that portion of the tree remains as it is
			
		step++;
	}
	
	cout << "Tree reconstruction in " << (int)step << " steps." << endl;
	
	//Print tree in file OutNew
	Print_tree_newick_form(phylo_tree,phylo_tree.begin(),OutNew);
	OutNew << ";\n";
	
    #if PRINT_TABLE
	//Close output file
	fclose(OutP);
    #endif

	//Close files
	fclose(InLCP);
    fclose(InCDA);
    fclose(InCDA_2);
    fclose(InEBWT);
	out_fb.close();
	
    fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	return 0;
}

