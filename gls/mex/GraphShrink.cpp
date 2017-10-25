#include "mex.h"
#include <cmath>
#include "blas.h"
//#include <matrix.h>

#define MAX_EDGE_NUM 400

typedef ptrdiff_t intt;

//void parse_params(const mxArray *params, intt *max_it, bool *pos);

 
typedef struct Edge_t
{   int id;
    int rmid;
    struct ELink_t * NodeLkArr[2];
} Edge;


 
typedef struct ELink_t
{
    struct Edge_t * edge;
    struct Node_t * node;
    struct ELink_t * next;
    struct ELink_t * prev;
} ELink;



typedef struct Node_t
{
    int    id;
    struct Data_t * DataList;
    struct ELink_t * EdgeList;
    struct Node_t * next;
    struct Node_t * prev;
} Node;


typedef struct Data_t
{
    int  Xid;
    struct Data_t * next;
} Data;


// graph add edge 
Edge * AddEdge(Node * Graph, Node * hnode, Node * tnode, int id)
{
    Edge * tedge = new Edge;
    ELink *  helk = new ELink;
    ELink *  telk = new ELink;
    
    tedge -> NodeLkArr[0] = helk;
    tedge -> NodeLkArr[1] = telk;
    tedge -> id = id;

    // helk
    helk -> edge = tedge;
    helk -> node = hnode;
    helk -> next = NULL; 
    
    ELink * elktemp = hnode -> EdgeList;
    
    if(NULL == elktemp)
    {
        hnode -> EdgeList = helk;
        helk -> prev = NULL;
    }
    else
    {
        while(NULL != elktemp->next)
        {
            elktemp = elktemp -> next;
        }
        elktemp -> next = helk;
        helk -> prev = elktemp;
    }
    
    //telk
    telk -> edge = tedge;
    telk -> node = tnode;
    telk -> next = NULL;
    
    elktemp = tnode -> EdgeList;
    if(NULL == elktemp)
    {
        tnode -> EdgeList = telk;
        telk -> prev = NULL;
    }
    else
    {
        while(NULL != elktemp->next)
        {
            elktemp = elktemp -> next;
        }
        elktemp -> next = telk;
        telk -> prev = elktemp;
    }
    
    return tedge;
}


//================ construct the Graph ===============
// node_list: xi  edge_id   edge_list: node_id, node_weight
Node * CreateGraph(int P, int nlM, double * node_list, int eM, double * edge_list, double * rmLab, Edge ** rmEdgeList)
{
    Node * Graph = NULL;
   
    Node ** NodeArr = new Node *[P];
    Node * nptr = new Node;
    nptr->EdgeList = NULL;
    nptr->next = NULL;
    nptr->prev = NULL;
    nptr->DataList = new Data;
    nptr->DataList->Xid = 1;
    nptr->DataList->next = NULL;
    nptr-> id = 1;
    NodeArr[0] = nptr;
    Graph = nptr;
    
    
    Node * preNode = nptr;
    // create node list
    for(int i = 1; i <P; i++)
    {
        nptr = new Node;
        nptr->EdgeList = NULL;
        nptr->next = NULL;
        nptr->prev = preNode;
        nptr->DataList = new Data;
        nptr->DataList->Xid = i+1;
        nptr->DataList->next = NULL;
        nptr-> id = i+1;
        preNode->next = nptr;
        preNode = nptr;
        NodeArr[i] = nptr;
    }

    // add edges
    int eBuf[3];
    Edge * tedge;
    int rmCnt = 0;
    for(int i = 0; i<eM; i++)
    {
        eBuf[0] = edge_list[i];
        eBuf[1] = edge_list[i + eM];
        eBuf[2] = edge_list[i + 2*eM];
       // eBuf[3] = edge_list[i + 3*eM];
        
        int head, tail;
        if(eBuf[2] > 0)
        {
            head = eBuf[0];
            tail =  eBuf[1];
        }
        else
        {
            head = eBuf[1];
            tail =  eBuf[0];
        }
        
       tedge =  AddEdge(Graph, NodeArr[head -1], NodeArr[tail -1], i+1);
       if(rmLab[i] == 1)
       {
          rmEdgeList[rmCnt] = tedge;
          tedge->rmid = rmCnt;
          rmCnt++;
       }
       else
       {
         tedge->rmid = -1;
       }
    }
    
    delete [] NodeArr;
    return Graph;
}

// release 
void ReleaseGraph(Node * Graph)
{
    Node * nptr = Graph;
    Node * tnptr;
    Data * dptr, * tdptr;
    ELink * elkptr, * telkptr;
    Edge * tedge;
    
    
    while(NULL != nptr)
    {
        tnptr = nptr -> next;
        
        dptr = nptr -> DataList;
        while(NULL != dptr)
        {
            tdptr = dptr -> next;
            delete dptr;
            dptr = tdptr;
        }
        
        elkptr = nptr-> EdgeList;
        while(NULL != elkptr)
        {
          telkptr = elkptr-> next;
          
          //==== delete the edge ====
          tedge = elkptr-> edge;
          if(tedge != NULL)
          {
              tedge-> NodeLkArr[0]->edge = NULL;
              tedge-> NodeLkArr[1]->edge = NULL;
              delete tedge;
          }
          //==========================
                  
          delete elkptr;
          elkptr = telkptr;
        }
        
        delete nptr;
        nptr = tnptr;
    }
}


// release 
void PrintGraph(Node * Graph)
{
    Node * nptr = Graph;
    Node * tnptr;
    Data * dptr, * tdptr;
    ELink * elkptr, * telkptr;
    Edge * tedge;
    
    
    while(NULL != nptr)
    {
        mexPrintf("\n====== i=%d ======\n", nptr->id);
        
        tnptr = nptr -> next;
        dptr = nptr -> DataList;
        mexPrintf("Data:\n");
        while(NULL != dptr)
        {
            tdptr = dptr -> next;
            mexPrintf("%d ", dptr->Xid);
            dptr = tdptr;
        }
        
        elkptr = nptr-> EdgeList;
        mexPrintf("\nEdge:\n");
        while(NULL != elkptr)
        {
          telkptr = elkptr-> next;
          
          //==== print the edge ====
          tedge = elkptr-> edge;
          mexPrintf("%d-%d  ", tedge-> NodeLkArr[0]->node->id, tedge-> NodeLkArr[1]->node->id);
          //==========================
          elkptr = telkptr;
        }
        
        nptr = tnptr;
        mexPrintf("\n");
    }
}



// check the edge list so that there is no duplicate edge
int MergeEdgeList(int & newuN, int * IdBuf, Edge ** EPtrBuf, ELink *  head, ELink * endhnode, int hid,  Edge ** rmEdgeList)
{
    int Nrm = 0;
    ELink *atmpelk, *btmpelk, *ctmpelk, *dtmpelk;
    int idbf = 0;
    int  tmpid;
            
    if(NULL != endhnode)
    {
        if(NULL  != endhnode->next)
        {
             atmpelk = head; // head list part
             while(atmpelk != endhnode->next)
                {
                   EPtrBuf[idbf] = atmpelk->edge;
                   tmpid = atmpelk->edge->NodeLkArr[0]->node->id;
                   if(tmpid != hid)
                   {
                      IdBuf[idbf] = tmpid;
                   }
                   else
                   {
                      IdBuf[idbf] = atmpelk->edge->NodeLkArr[1]->node->id;
                   }      
                   idbf++;
                   atmpelk = atmpelk-> next;
                }
             
              btmpelk = endhnode;  // tail list part
              while(NULL != atmpelk)
               {
                   int loct = 0;
                   tmpid = atmpelk->edge->NodeLkArr[0]->node->id;
                   if(tmpid == hid)
                   {
                     tmpid = atmpelk->edge->NodeLkArr[1]->node->id;
                     loct = 1;
                   }    
                   
                   int flag = -1;
                   for(int ii = 0; ii < idbf; ii++)
                   {
                       if(tmpid == IdBuf[ii])
                       {
                           flag = ii;
                       }
                   }
                   
                   // merge the duplicate edges 
                   if(flag > -1)
                   {
                       btmpelk->next = atmpelk->next;
                       // skip the node
                       if(NULL != atmpelk->next)
                       {
                          atmpelk->next->prev = btmpelk;              
                       } 
                       
                       ctmpelk = atmpelk;
                       atmpelk = atmpelk->next;
                       dtmpelk = ctmpelk->edge->NodeLkArr[loct];
                       
                       if(NULL != dtmpelk->next)
                       {
                              dtmpelk->next->prev = dtmpelk->prev;
                        }
                       
                       if(NULL == dtmpelk->prev)
                       {
                           dtmpelk->node->EdgeList = dtmpelk->next;
                       }
                       else
                       {
                           dtmpelk->prev->next = dtmpelk-> next;
                       }
                       
                       if( ctmpelk->edge-> rmid >= 0)
                       {
                           if( EPtrBuf[flag]-> rmid <0 )
                           {
                             EPtrBuf[flag]-> rmid  = ctmpelk->edge-> rmid;
                             rmEdgeList[ctmpelk->edge-> rmid] = EPtrBuf[flag];
                           }
                           else
                           {
                             rmEdgeList[ctmpelk->edge-> rmid] = NULL;
                           }
                       }
                       
                       
                     //  EdgePtrList[ctmpelk->edge->id -1] = EPtrBuf[flag];
                       delete dtmpelk;
                       delete ctmpelk->edge;
                       delete ctmpelk;
                       Nrm++;
                       newuN --;
                   }
                   else
                   {
                     btmpelk = btmpelk -> next;
                     atmpelk = atmpelk -> next;
                   }
               }
        }
    }

    return Nrm;
}

//rmM, rmLab, EdgePtrList
void RemoveEdge(int & newP, int & newuN, Node *& G, int rmN,  Edge ** rmEdgeList)
{
   Edge * tedge;
   Node * hnode, * tnode;
   ELink * helk, * telk, *npelk, *endhnode;
   Data * datptr;
 
   int *IdBuf = new int[MAX_EDGE_NUM]; 
   Edge ** EPtrBuf = new Edge *[MAX_EDGE_NUM];
   int hid, tid, tmpid;
   int rmEId;
   
   
   
    for(int i = 0; i<rmN; i++)
    {

        tedge = rmEdgeList[i];
        
        if(NULL == tedge)
        {
          continue;
        }
        
          mexPrintf("entered for loop i=%d\n", i);
        
   //     mexPrintf("\n rm_edge = %d\n", i);
 
        helk = tedge-> NodeLkArr[0];
        mexPrintf("step01\n");
        telk = tedge-> NodeLkArr[1];
          mexPrintf("step02\n");
          
          if(NULL == helk)
          {
               mexPrintf("helk is null\n");
          }
          
            if(NULL == telk)
          {
               mexPrintf("helk is null\n");
          }
          
          
        hnode = helk -> node;
          mexPrintf("step03\n");
        tnode = telk -> node;
          mexPrintf("step04\n");
        
        hid = hnode-> id;
             mexPrintf("step05\n");
        tid = tnode-> id;
        
        mexPrintf("step1\n");
        
    //    mexPrintf("\n ===== hid=%d, tid=%d ======\n", hid, tid);
        
        ///=========== merge the data list ========== //
        datptr = hnode-> DataList;
        while(datptr->next != NULL)
        {
            datptr = datptr->next;
        }
        datptr->next = tnode ->DataList;
        
        mexPrintf("step2\n");
        
        //======= remove helk ====================///
        if(helk->prev == NULL)
        {
            hnode->EdgeList = helk->next;
        }
        else
        {
            helk->prev->next = helk->next;
        }
        
        if(NULL !=  helk->next)
        {
            helk->next->prev = helk->prev;
        }
        
        mexPrintf("step3\n");
       
        //===== remove telk =======================//
        if(telk -> prev == NULL)
        {
            tnode->EdgeList = telk->next;
        }
        else
        {
            telk -> prev -> next = telk-> next;
        }
        
        // ptr for previous
        if(NULL != telk->next)
        {
            telk->next->prev = telk->prev;
        }
        
        mexPrintf("step4\n");
        //==================merge the edge lists ===========//
        endhnode = hnode->EdgeList;
        if(NULL == endhnode)
        {
           hnode->EdgeList = tnode->EdgeList;
        }
        else
        {
            while(NULL !=   endhnode->next)
            {
                 endhnode = endhnode -> next;
            }
            endhnode -> next = tnode->EdgeList;
            if(NULL != tnode->EdgeList)
            {
              tnode->EdgeList->prev = endhnode;
            }
        }
        
        mexPrintf("step5\n");
    //    mexPrintf("\n done, going to merge !\n");
        
        if(NULL != endhnode)
        {
            npelk =  endhnode -> next; 
            while(NULL != npelk)
            {
                npelk->node = hnode;
                npelk = npelk->next;
            }
        }
        
        
//         ELink * tmpptr = hnode->EdgeList;
//         while(NULL != tmpptr)
//         {
//             mexPrintf(" %d-%d ", tmpptr->edge->NodeLkArr[0]->node->id, tmpptr->edge->NodeLkArr[1]->node->id);
//             tmpptr = tmpptr->next;
//         }
        
          mexPrintf("step6\n");
        
       MergeEdgeList(newuN, IdBuf, EPtrBuf, hnode->EdgeList,endhnode, hid, rmEdgeList);

//         mexPrintf(" \n");
//         tmpptr = hnode->EdgeList;
//         while(NULL != tmpptr)
//         {
//             mexPrintf(" %d-%d ", tmpptr->edge->NodeLkArr[0]->node->id, tmpptr->edge->NodeLkArr[1]->node->id);
//             tmpptr = tmpptr->next;
//         }
         mexPrintf("step7\n");
        //========================== delete the tnode=============================== //
       if(NULL == tnode->prev)
       {
            G = tnode->next;
       }
       else
       {
            tnode->prev->next = tnode->next;   
       }
      
           mexPrintf("step8\n");
         
        if(tnode->next != NULL)
        {
            tnode->next-> prev = tnode->prev;
        }
        
   //     mexPrintf("\n node_id = %d\n", tnode->id);
     //   EdgePtrList[tedge->id-1] = NULL; 
        delete tedge;
        delete helk;
        delete telk;
        delete tnode;
        
        newP--;
        newuN--;
        
     //   mexPrintf("\n ^^^^^^^^^^===================^^^^^^^^ \n", tnode->id);
     //   PrintGraph(G);
    }
   
    delete [] IdBuf;
    delete [] EPtrBuf;
}


///////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

// Input: X,D,A,lambda, max_iter, optD
// bool debug = true;
//  if(debug)
//  {
//     mexPrintf("Entered!\n");
//  }
 
 intt rmM = mxGetM(prhs[0]); 
 intt nlM = mxGetM(prhs[1]);
 intt nlN = mxGetN(prhs[1]);
 intt eM =  mxGetM(prhs[2]);
 intt eN =  mxGetN(prhs[2]);
 
 double *rmLab, *node_list, *edge_list; //*Q, *lambda, *alpha, *max_iter, *OptThr, *D, *SigmaSqu;

 rmLab = (double*)mxGetPr(prhs[0]);
 node_list = (double*)mxGetPr(prhs[1]);
 edge_list = (double*)mxGetPr(prhs[2]);
 int P = node_list[nlM-1];
 
//   if(debug)
//  {
//     mexPrintf("rmM=%d, nlM=%d,  nlN=%d, eM=%d,  eN=%d, P = %d\n", rmM, nlM, nlN, eM, eN, P);
//   }
 
 int rmN = 0;
 for(int i=0; i<rmM; i++)
 {
     if(rmLab[i] == 1)
     {
         rmN ++;
     }
 }
 
 Edge ** rmEdgeList = new Edge *[rmN];
 Node* G = CreateGraph(P, nlM, node_list, eM, edge_list, rmLab, rmEdgeList);
 //PrintGraph(G);
// mexPrintf("=================XXXXXXXX=======================");
 int newP = P;
 int newuN = eM;
 
 RemoveEdge(newP, newuN, G, rmN, rmEdgeList);
 
 plhs[0]  = mxCreateDoubleMatrix(P, 2, mxREAL);
 double * OutNode = (double *)mxGetPr(plhs[0]);
 plhs[1]  = mxCreateDoubleMatrix(newuN, newP, mxREAL);
 double * newD = (double *)mxGetPr(plhs[1]);
 
// PrintGraph(G);
 
 Node * nptr = G;
 Data * dptr;
 ELink * elkptr;
 Edge * eptr;
         
 //int nid=1;
 int i = 0;
 int nid = 1;
 
 while(NULL != nptr)
 {
     nptr->id = nid;
     nptr = nptr -> next;
     nid++;
 }
 
 nptr = G;
 //nid = 1;
 int eid = 0;
 int id1, id2;
 while(NULL != nptr)
 {
     dptr = nptr->DataList;
     elkptr = nptr->EdgeList;
     
     while(NULL!= dptr)
     {
         OutNode[i] = nptr->id ;
         OutNode[i+P] = dptr->Xid;
         dptr = dptr->next;
         i++;
     }
     
     while(NULL != elkptr)
     {
        eptr = elkptr->edge;
        
        if(-2 != eptr-> rmid)
        {
            id1 = eptr->NodeLkArr[0]->node->id;
            id2 = eptr->NodeLkArr[1]->node->id;

            newD[eid + newuN*(id1-1)] = 1; 
            newD[eid + newuN*(id2-1)] = -1; 
            eptr-> rmid = -2;
            eid++;
       //      mexPrintf("nid=%d, %d-%d\n", nptr->id, id1,id2);
        }
                
        elkptr = elkptr->next;
     }

     nptr = nptr -> next;
   //  nid ++;
 }

 //PrintGraph(G);
 
 ReleaseGraph(G);
 delete [] rmEdgeList;

}



