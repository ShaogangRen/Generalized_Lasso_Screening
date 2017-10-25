function T = get_transformationT(P, Node, Edge)


ANode = cell(P,1);
%SNode = ones(P,1);
FNodeE = cell(P,1);
VFlag = zeros(P,1);

for i = 1:P
    VFlag(i) = length(Node{i});
    ANode{i} = i;
    FNodeE{i} = true(1,length(Node{i}));
end


roundN = 0;
while 1
    roundN = roundN +1;
    sthupdated = false;
  %  disp(['====================round=' num2str(roundN) '======================' ]);
    
    for i = 1:P
      %  disp(['$$$$$$$$$$$$ i=' num2str(i) '$$$$$$$$$$$$$$' ]);
        if VFlag(i) == 1
            sthupdated = true;
            
            nbuf = [ ];
            run = true;
            ni = i;
            fromE = 0;

            while run
                ANode{ni} = [ANode{ni}, nbuf];
                nbuf = ANode{ni};
                VFlag(ni) = VFlag(ni)  -1;
                if fromE > 0
                   FNodeE{ni}(Node{ni} == fromE) = false;  
                end

                if VFlag(ni) > 1
                   run = false;
                elseif any(FNodeE{ni})
                    
                     EId = Node{ni}(FNodeE{ni});
                     FNodeE{ni}(Node{ni} == EId) = false;

                     Edg =  Edge(EId,:);
                     ni = Edg(Edg ~= ni); %% new node id
                     fromE = EId;   
                else
                     run = false;
                end
            end
        end
    end
    
    if ~ sthupdated
        break;
    end
end

T = [];
for i = 1:P
    tv = zeros(P,1);
    tv(ANode{i}) = 1;
    T = [T, tv];
end





