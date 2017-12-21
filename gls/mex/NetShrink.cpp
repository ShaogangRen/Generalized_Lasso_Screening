/* By Shaogang Ren, renshaogang@gmail.com. All rights reserved. */
#include "mex.h"
#include <cmath>
#include "blas.h"
//#include <matrix.h>

#define MAX_EDGE_NUM 200

typedef ptrdiff_t intt;

typedef struct Edge_t
{   int id;
    int rm;//1: remove, 0: not remove
    int state; // 1: active, 0: inactive
    int Nodes[2];
} Edge;

 
typedef struct ELink_t
{
    int edge;
    struct ELink_t * next;
} ELink;


typedef struct Node_t
{
    int    id;
    int    state; // 0: nonactive, 1: active
    struct Data_t * DataList;
    struct ELink_t * EdgeList;
} Node;


typedef struct Data_t
{
    int  Xid;
    struct Data_t * next;
} Data;


void BuildEdgeList(int eM, double * edge_list, Edge * EdgeList, Node * NodeList, double * rmLab, int *rmEdgeList)
{
     int eBuf[3];
     Node * tnode;
     ELink * telk, *mvelk;
     int rmidx = 0;
     
     for(int i=0; i<eM; i++)
     {
       eBuf[0] = edge_list[i];
       eBuf[1] = edge_list[i + eM];
       eBuf[2] = edge_list[i + 2*eM];        

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

      //  mexPrintf("\n add edgei=%d, head=%d, tail=%d \n",i, head, tail);
       EdgeList[i].id = i;  
       EdgeList[i].state = 1;
       EdgeList[i].Nodes[0] = head -1;
       EdgeList[i].Nodes[1] = tail -1;

       // mexPrintf("\n add edgei=%d, head=%d, tail=%d \n",i, EdgeList[i].Nodes[0], EdgeList[i].Nodes[1]);
       
       if(1== rmLab[i])
       {
           rmEdgeList[rmidx] = i;
           EdgeList[i].rm = rmidx;
           rmidx++;
       }
       else
       {
           EdgeList[i].rm = -1;
       }
       
       tnode = &NodeList[EdgeList[i].Nodes[0]];
       telk = new ELink;
       telk->edge = i;
       telk->next = NULL;
       
       mvelk = tnode->EdgeList;
       if(NULL == mvelk)
       {
         tnode->EdgeList = telk;
       }
       else
       {
          while(NULL != mvelk->next)
          {
              mvelk = mvelk->next;
          }
          mvelk->next = telk;
       }
       
       tnode = &NodeList[EdgeList[i].Nodes[1]];
       telk = new ELink;
       telk->edge = i;
       telk->next = NULL;
       
       mvelk = tnode->EdgeList;
       if(NULL == mvelk)
       {
         tnode->EdgeList = telk;
       }
       else
       {
          while(NULL != mvelk->next)
          {
              mvelk = mvelk->next;
          }
          mvelk->next = telk;
       }
       
     }
}


void ReleaseGraph(int P, Node* NodeList, Edge * EdgeList)
{
    Data * dptr, * tdptr;
    ELink * elkptr, * telkptr;
    
    for(int i=0; i<P; i++)
    {
       if(1 == NodeList[i].state)
       {
            dptr = NodeList[i].DataList;
            
            while(NULL != dptr)
            {
              tdptr = dptr->next;
              delete dptr;
              dptr = tdptr;
            }
            
            elkptr = NodeList[i].EdgeList;
            while(NULL != elkptr)
            {
               telkptr = elkptr->next;
               delete elkptr;
               elkptr = telkptr; 
            }
       } 
    }
    
    delete [] NodeList;
    delete [] EdgeList;
}

void PrintGraph(int P, Node* NodeList, Edge * EdgeList)
{
    Data * dptr, * tdptr;
    ELink * elkptr;
    
    for(int i=0; i<P; i++)
    {
       if(1 == NodeList[i].state)
       {
            mexPrintf("\n====== i=%d ======\n", (NodeList[i].id));
            dptr = NodeList[i].DataList;
            
            mexPrintf("Xi:\n");
            while(NULL != dptr)
            {
              mexPrintf("%d ", dptr->Xid);
              dptr = dptr->next;
            }
            
            mexPrintf("\nEdge: \n");
            elkptr = NodeList[i].EdgeList;
            while(NULL != elkptr)
            {
               mexPrintf("%d-%d ", NodeList[EdgeList[elkptr->edge].Nodes[0]].id, NodeList[EdgeList[elkptr->edge].Nodes[1]].id);
               elkptr = elkptr->next;
            }
       }   
    }
    
    mexPrintf("\n");
}

ELink * MergeEdgeList(Node* NodeList, Edge * EdgeList, int *rmEdgeList, ELink * hlist, ELink * tlist, int hnode, int tnode)
{
   ELink * elkptr, * telkptr;
   int tonode;
   ELink * resptr, *aptr, *bptr, *cptr, *dptr;
   int idxt =0;
   int nid1, nid2;
   int tnidx, rnidx;
   
   resptr = NULL;
        
     elkptr = hlist;
     while(NULL != elkptr)
     {
        nid1 =  EdgeList[elkptr->edge].Nodes[0];
        nid2 =  EdgeList[elkptr->edge].Nodes[1];
        tonode = ((nid1 == hnode) ? nid2 : nid1);
     
        if(tonode != tnode)
        {
//             IdBuf[idxt]  = tonode; 
//             EdgeBuf[idxt] = elkptr->edge;
//             idxt ++;
            
            if(NULL == resptr)
            {
               resptr = elkptr;
               aptr = elkptr;
            }
            else
            {
               aptr->next = elkptr;
               aptr = elkptr;
            }
            
            elkptr = elkptr->next;
        }
        else
        {
            telkptr = elkptr;
            elkptr = elkptr->next;
            EdgeList[telkptr->edge].state = 0;
            delete telkptr;
        }
     }
     
     // the tail list
     elkptr = tlist;
     while(NULL != elkptr)
     {
        nid1 =  EdgeList[elkptr->edge].Nodes[0];
        if(tnode == nid1)
        {
            tnidx = 0; tonode = EdgeList[elkptr->edge].Nodes[1];
        }
        else
        {
            tnidx = 1; tonode = nid1;
        }
        
        bool flag = true;
    //    int  pidx = -1;
        
        //determine if can be reused
        if(tonode == hnode)
        {
          flag = false;
        }

        if(flag)
        {
                EdgeList[elkptr->edge].Nodes[tnidx] = hnode;
            
                if(NULL == resptr)
                {               
                   resptr = elkptr;
                   aptr = elkptr;
                }
                else
                {
                   aptr->next = elkptr;
                   aptr = elkptr;
                }
                elkptr = elkptr->next;
        }
        else
        { 
            telkptr = elkptr;
            elkptr = elkptr->next;
            
            EdgeList[telkptr->edge].state = 0;
            delete telkptr;
        }
     }
     
     aptr->next = NULL;
     return resptr;
}


int RemoveEdge(int P, Node* NodeList, Edge * EdgeList, int rmN, int * rmEdgeList)
{
    int rmid;
    Edge * rmEg;
    int hnode, tnode;
    Data * dptr;
    
   // int *IdBuf = new int[MAX_EDGE_NUM]; 
   // int *EdgeBuf = new int[MAX_EDGE_NUM]; 
    int ididx;
    ELink * elkptr;
    int nid1, nid2;
    int newP = P;
    
   //  int tnode, hnode;
 
    for(int i=0; i<rmN; i++)
    {
         rmid = rmEdgeList[i];
        //EdgeList[rmid];
         
        if(1 == EdgeList[rmid].state)
        {
             hnode = EdgeList[rmid].Nodes[0];
             tnode = EdgeList[rmid].Nodes[1];
        //     mexPrintf(" remove edge: %d: %d-%d \n", rmid, hnode, tnode); 
             dptr = NodeList[hnode].DataList;
             while(NULL != dptr->next)
             {
                 dptr = dptr->next;
             }
             dptr->next = NodeList[tnode].DataList;
             
             NodeList[hnode].EdgeList =  MergeEdgeList(NodeList, EdgeList, rmEdgeList, NodeList[hnode].EdgeList, NodeList[tnode].EdgeList, hnode, tnode);

             NodeList[tnode].state = 0;
             newP --;
             EdgeList[rmid].state = 0;
        }
    }
    
   // delete [] IdBuf;
   // delete [] EdgeBuf;
    return newP;
}

///////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

 intt rmM = mxGetM(prhs[0]); 
 intt nlM = mxGetM(prhs[1]);
 intt nlN = mxGetN(prhs[1]);
 intt eM =  mxGetM(prhs[2]);
 intt eN =  mxGetN(prhs[2]);
 
 double *rmLab, *node_list, *edge_list; 

 rmLab = (double*)mxGetPr(prhs[0]);
 node_list = (double*)mxGetPr(prhs[1]);
 edge_list = (double*)mxGetPr(prhs[2]);
 int P = node_list[nlM-1];
 int newP;

 Node * NodeList = new Node[P];
 Edge * EdgeList = new Edge[eM];

 
 // remove edge list
 //mexPrintf("\n rmM= %d \n", rmM);
 int rmcnt = 0;
 for(int i=0; i<rmM; i++)
 {
     if(1 == rmLab[i])
     {
         rmcnt ++;     
     }
 }
 
 // construct the node list
 for(int i= 0; i< P; i++)
 {
     NodeList[i].id = i;
     NodeList[i].state = 1;
     NodeList[i].DataList = new Data;
     NodeList[i].DataList->Xid = i+1;
     NodeList[i].DataList->next = NULL;
     NodeList[i].EdgeList = NULL;
 }

//construct the edge list
 int *rmEdgeList = new int[rmcnt];
 
 BuildEdgeList(eM, edge_list, EdgeList, NodeList, rmLab, rmEdgeList);
 
 newP = RemoveEdge(P, NodeList, EdgeList, rmcnt, rmEdgeList);
 
 //mexPrintf("\n === rmcnt= %d \n", rmcnt);
 //mexPrintf("\n =========================================================== \n");        
// mexPrintf("\n newP= %d \n", newP);
 
 // PrintGraph(P, NodeList, EdgeList);
   
 plhs[0]  = mxCreateDoubleMatrix(P, 2, mxREAL);
 double * OutNode = (double *)mxGetPr(plhs[0]);
 int * A = new int[newP*newP];
 for(int i=0;i<newP;i++)
    for(int j=0; j<newP;j++)
        A[i + j*newP] = 0;
 
 int nidx = 0;
 int idx = 0;
 int newuN = 0;
 int id1, id2;
 Data * dptr;
 ELink * elkptr;
 
 for(int i=0; i<P; i++)
 {
    if(1 == NodeList[i].state)
    {
        NodeList[i].id = nidx++;
    }
 }
 
//  PrintGraph(P, NodeList, EdgeList);
 
 //mexPrintf("\n nidx=%d, newP= %d \n", nidx, newP);
 
 for(int i=0; i<P; i++)
 {
       if(1 == NodeList[i].state)
       {
          dptr = NodeList[i].DataList;
          elkptr =  NodeList[i].EdgeList;

          while(NULL!= dptr)
          {
            OutNode[idx] = NodeList[i].id+1;
            OutNode[idx + P] =  dptr->Xid;
            idx ++;
            dptr = dptr->next;
          }
        
          while(NULL != elkptr)
          {
            if(EdgeList[elkptr->edge].rm != -2)
            {
                EdgeList[elkptr->edge].rm = -2;
                id1 = NodeList[EdgeList[elkptr->edge].Nodes[0]].id;
                id2 = NodeList[EdgeList[elkptr->edge].Nodes[1]].id;
                
                if(0 == A[id1+id2*newP])
                {
                      newuN++;
                }
                
                A[id1+id2*newP]++; A[id2+id1*newP]++;
          
                
            }
            elkptr = elkptr->next;
          }
       }
 }
 
// mexPrintf("\n newuN = %d \n", newuN); 
 
plhs[1]  = mxCreateDoubleMatrix(newuN, newP, mxREAL);
double * newD = (double *)mxGetPr(plhs[1]);
 
int ui = 0;
for(int i= 0; i<newP; i++)
   for(int j= i+1; j< newP; j++)
   {
        if(A[i+j*newP] > 0)
        {
          newD[ui + i*newuN] = A[i+j*newP];
          newD[ui + j*newuN] = - A[i+j*newP];
          ui++;
        }
   }
 
 ReleaseGraph(P, NodeList, EdgeList);
 delete [] rmEdgeList;
 delete [] A;
}






