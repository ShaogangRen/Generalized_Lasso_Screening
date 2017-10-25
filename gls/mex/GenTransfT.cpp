#include "mex.h"
#include <cmath>
#include "blas.h"
//#include <matrix.h>

//#define MAX_EDGE_NUM 200

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
    int visited;
    struct ELink_t * next;
} ELink;


typedef struct Node_t
{
    int    id;
    int    state; // 0: nonactive, 1: active
    int    ActEgNum;
    int    PathNum;
    struct Data_t * DataList;
    struct ELink_t * EdgeList;
} Node;


typedef struct Data_t
{
    int  Xid;
    struct Data_t * next;
} Data;

typedef struct LNode_t
{
 struct Data_t * DataList;
 struct LNode_t * next; 
} LNode;



void BuildEdgeList(int eM, double * edge_list, Edge * EdgeList, Node * NodeList)
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

       // head node
       tnode = &NodeList[EdgeList[i].Nodes[0]];
       telk = new ELink;
       telk->edge = i;
       telk->visited = 0;
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
       tnode->ActEgNum ++;
       
       // tail node
       tnode = &NodeList[EdgeList[i].Nodes[1]];
       telk = new ELink;
       telk->edge = i;
       telk->visited = 0;
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
       tnode->ActEgNum ++;
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
            //mexPrintf("\n====== i=%d ======\n", (NodeList[i].id));
            dptr = NodeList[i].DataList;
            
           // mexPrintf("Xi:\n");
            while(NULL != dptr)
            {
             // mexPrintf("%d ", dptr->Xid);
              dptr = dptr->next;
            }
            
           // mexPrintf("\n ActEgNum=%d, PathNum=%d ", NodeList[i].ActEgNum, NodeList[i].PathNum);

          //  mexPrintf("\nEdge: \n");
            elkptr = NodeList[i].EdgeList;
            while(NULL != elkptr)
            {
             //  mexPrintf("%d-%d ", NodeList[EdgeList[elkptr->edge].Nodes[0]].id, NodeList[EdgeList[elkptr->edge].Nodes[1]].id);
               elkptr = elkptr->next;
            }
       }   
    }
    
   // mexPrintf("\n");
}


void ConstrSgNdList(int P, int eM,Edge * EdgeList, Node* NodeList)
{

    bool sthupdated = false;
    bool run;
    Data * dhead,* dtail, *dptr, *tdptr;
    ELink * telk;
    Node * tNdPtr;
    int tnd0, from, nextid;

    while(true)
    {
        sthupdated = false;
        for(int i= 0; i<P; i++)
        {
            if(1 == NodeList[i].ActEgNum)
            {
                sthupdated = true;
                run = true;
                
                dhead = NULL;
                dtail = NULL;
                tNdPtr = &NodeList[i];
                
                from = -1;
                
                while(run)
                {
                  if(NULL != dhead)
                  {
                       dtail ->next = tNdPtr->DataList ;
                       tNdPtr->DataList = dhead;
                  }
                 
                  // copy the datalist
                   dhead = new Data;
                   dhead ->Xid = tNdPtr->DataList->Xid;
                   dtail = dhead;
                   dptr = tNdPtr->DataList->next;
                   while(NULL != dptr)
                   {
                         tdptr = new Data;
                         tdptr->Xid = dptr->Xid;
                         tdptr->next = NULL;
                         dtail -> next = tdptr;
                         dtail = tdptr;
                         dptr = dptr->next;
                   }
                 
                  // set the node state 
                   if(from >-1)
                   {
                      telk = tNdPtr->EdgeList;
                      while(EdgeList[telk->edge].id != from)
                      {
                        telk = telk->next;
                      }
                      telk->visited = 1;
                   }
                   
                   
                  tNdPtr->ActEgNum --;
                  if(tNdPtr->ActEgNum > 1)
                  {
                     run = false;
                  }
                  else 
                  {
                      telk = tNdPtr->EdgeList;
                      bool conti = false;
                      
                      while(NULL != telk)
                      {
                        if( telk->visited == 0)
                        {
                            telk->visited  = 1;
                            from =  telk->edge;
                            tnd0 = EdgeList[telk->edge].Nodes[0];
                            if(tnd0 == tNdPtr->id)
                            {
                              nextid = EdgeList[telk->edge].Nodes[1];
                            }
                            else
                            {
                              nextid = tnd0;   
                            }
        
                            tNdPtr = &NodeList[nextid];
                            conti = true;
                        }
                        telk = telk->next;
                      }
                      run = conti;
                  }
                }
                
            }
        }
        
        
        if(!sthupdated)
        {
         break;   
        }
    }
}

LNode * TriNdLink(int P, int eM, Node * NodeList, Edge * EdgeList, int & ndNum)
{
   //for()
    int trnd;
    // choose one node with the number of active edge >3
    for(int i = 0; i<P; i++)
    {
        if(NodeList[i].ActEgNum > 2)
        {
          trnd = i;
          break;   
        }
    }
    
    LNode * ListTable = new LNode;
    ListTable->DataList = NULL;
    ListTable->next = NULL;
    
    LNode * preNd = ListTable;
    ELink * telk, *nelk; 
    int tnd0, nextid, fromEg;
    bool run = true;
    
    ndNum = 0;
    
    while(run)
    {
        
        if(NodeList[trnd].PathNum == 0)
            break;   
        
        LNode * tlnode = new LNode;
        preNd->next = tlnode;
        tlnode->next = NULL;
        
        Data * dptr = new Data;
        dptr->Xid = trnd;
        ndNum++;
        dptr -> next = NULL;
        tlnode->DataList = dptr;
        Data * predptr = dptr;
        

        telk = NodeList[trnd].EdgeList;
        while(NULL != telk)
        {
            if(telk->visited == 0)
            {
                break;
            }
            telk = telk->next;
        }
        telk->visited = 1;
        NodeList[trnd].PathNum --;

        ////////////////
        tnd0 = EdgeList[telk->edge].Nodes[0];
        if(tnd0 == trnd)
        {
          nextid = EdgeList[telk->edge].Nodes[1];
        }
        else
        {
          nextid = tnd0;   
        }
        trnd = nextid;
        fromEg = telk->edge;
        
        while( NodeList[trnd].ActEgNum == 2)
        {
            telk = NodeList[trnd].EdgeList;
            while(NULL != telk)
            {
                if(telk-> edge == fromEg)
                {
                    break;
                }
                telk = telk->next;
            }
            telk->visited = 1;
            NodeList[trnd].PathNum --;
             
           dptr = new Data;
           dptr->Xid = trnd;
           ndNum++;
           dptr -> next = NULL;
           predptr->next = dptr;
           predptr = dptr;

            telk = NodeList[trnd].EdgeList;
            while(NULL != telk)
            {
                if(telk->visited == 0)
                {
                    break;
                }
                telk = telk->next;
            }
            telk->visited = 1;
            NodeList[trnd].PathNum --;

            ////////////////
            tnd0 = EdgeList[telk->edge].Nodes[0];
            if(tnd0 == trnd)
            {
              nextid = EdgeList[telk->edge].Nodes[1];
            }
            else
            {
              nextid = tnd0;   
            }
            trnd = nextid;
            fromEg = telk->edge;
         }

        telk = NodeList[trnd].EdgeList;
        while(NULL != telk)
        {
            if(telk-> edge == fromEg)
            {
                break;
            }
            telk = telk->next;
        }
        telk->visited = 1;
        NodeList[trnd].PathNum --;
          
       
       dptr = new Data;
       dptr->Xid = trnd;
       ndNum++;
       dptr -> next = NULL;
       predptr->next = dptr;
      // predptr = dptr;
        
      preNd = preNd->next;
    }
    
    return ListTable;
}


void ReleaseTrNdList(LNode * LNList)
{
    LNode * tlnode, *nlnode;
    Data * tdata, *ndata;
    tlnode = LNList->next;
    delete LNList;
    
    while(NULL != tlnode)
    {
        tdata = tlnode->DataList;
        while(NULL != tdata)
        {
            ndata = tdata->next;
            delete tdata;
            tdata = ndata;
        }
        
       nlnode = tlnode->next;
       delete tlnode;
       tlnode = nlnode;
    }  
}

void PrintTrNdList(LNode * LNList)
{
    LNode * tlnode;
    Data * tdata;
    tlnode = LNList->next;
   
    //mexPrintf("\n+++++++++++++++++++++++++\n");
    while(NULL != tlnode)
    {
        tdata = tlnode->DataList;
        while(NULL != tdata)
        {
            //mexPrintf("%d ", tdata->Xid);
            tdata = tdata->next;
        }
        //mexPrintf("\n");
       tlnode = tlnode->next;
    }  
}

int NodeDataCnt(int P, Node * NodeList)
{
    int ndcnt = 0;
    Data * tdata;
    for(int i=0; i<P;i++)
    {
        tdata = NodeList[i].DataList;
        while(NULL != tdata)
        {
            ndcnt++;
            tdata = tdata->next;
        }
    }
    
    return ndcnt;
}

///////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
 double  *NdNum, *edge_list; 

 intt eM =  mxGetM(prhs[1]);
 intt eN =  mxGetN(prhs[1]);

 NdNum = (double*)mxGetPr(prhs[0]);
 edge_list = (double*)mxGetPr(prhs[1]);
 int P = NdNum[0];
 int newP;

 Node * NodeList = new Node[P];
 Edge * EdgeList = new Edge[eM];
 
 // remove edge list
 mexPrintf("\n eM= %d, eN= %d, P= %d \n", eM, eN, P);

 // construct the node list
 for(int i= 0; i< P; i++)
 {
     NodeList[i].id = i;
     NodeList[i].state = 1;
     NodeList[i].DataList = new Data;
     NodeList[i].DataList->Xid = i+1;
     NodeList[i].DataList->next = NULL;
     NodeList[i].EdgeList = NULL;
     NodeList[i].ActEgNum = 0;
     NodeList[i].PathNum = 0;
 }

//construct the edge list
 BuildEdgeList(eM, edge_list, EdgeList, NodeList);
 PrintGraph(P, NodeList, EdgeList);
 ConstrSgNdList(P, eM,EdgeList, NodeList);
 
 int NdDataNum = NodeDataCnt(P, NodeList); 

 for(int i = 0; i<P; i++)
 {
    NodeList[i].PathNum = NodeList[i].ActEgNum;
 }
 
 int ListNdNum;
 LNode * ListNodeTable =  TriNdLink(P, eM, NodeList, EdgeList, ListNdNum);
 
 mexPrintf("\n NdDataNum=%d, ListNodeNum=%d \n", NdDataNum, ListNdNum );
 plhs[0]  = mxCreateDoubleMatrix(NdDataNum, 2, mxREAL);
 double * NodeData = (double *)mxGetPr(plhs[0]);

 plhs[1]  = mxCreateDoubleMatrix(ListNdNum, 2, mxREAL);
 double * ListNode = (double *)mxGetPr(plhs[1]);
 
 

Data * tdata;
int idxi =0;
for(int i=0; i<P;i++)
{
    tdata = NodeList[i].DataList;
    while(NULL != tdata)
    {
        NodeData[idxi] = i+1;
        NodeData[idxi + NdDataNum] = tdata->Xid;
        idxi ++;
        tdata = tdata->next;
    }
}
    

    LNode * tlnode;
    tlnode = ListNodeTable->next;
    idxi = 0; int listcnt = 0;
//    int st;
    while(NULL != tlnode)
    {
        tdata = tlnode->DataList;
//         if(NULL != tdata)
//         {
//             ListNode[idxi] = listcnt + 1;
//             ListNode[idxi + ListNdNum] = tdata->Xid + 1;
//             st = tdata->Xid + 1;
//             idxi++;
//             tdata = tdata->next; 
//         }
        
        while(NULL != tdata)
        {
            ListNode[idxi] = listcnt + 1;
            ListNode[idxi + ListNdNum] = tdata->Xid + 1;
            idxi++;
            tdata = tdata->next;
        }
        
       listcnt ++;
       tlnode = tlnode->next;
    }  
 
 
 
 //mexPrintf("\n============================================\n");
 
 PrintTrNdList(ListNodeTable);
 ReleaseTrNdList(ListNodeTable);
 
 PrintGraph(P, NodeList, EdgeList);
 ReleaseGraph(P, NodeList, EdgeList);

}






