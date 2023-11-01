clc; %Clear command window
clear all;
%t0 = clock; %Current date and time as date vector
mkdir dataWPN; %create fold storing data files


allorthology='Nucleus_protein_allorthology.txt';
DATAfile =allorthology;
fid = fopen(DATAfile,'r'); %open the file 
Df = textscan(fid,'%s%s%d%s'); %read file
fclose(fid); %close file, success 0, fail -1 
protein=cat(2,Df{1},Df{2},Df{4});
essentialproteinnumbers = sum(Df{3});
PIN='krogan2006_extended.txt';
DATAfile_r=PIN;
fid_r = fopen(DATAfile_r,'r'); %open file
Df_r = textscan(fid_r,'%s%s'); %read file
fclose(fid_r); %close file success 0,fail -1
Pro = union(Df_r{1},Df_r{2}); 
number = length(Pro); 
IN='yeastgenedata.txt';
fin = fopen(IN,'r'); %open file
Din = textscan(fin,'%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'); %read file
fclose(fin); %close file success 0,fail -1
express=cat(2,Din{3},Din{4},Din{5},Din{6},Din{7},Din{8},Din{9},Din{10},Din{11},Din{12},Din{13},Din{14},Din{15},Din{16},Din{17},Din{18},Din{19},Din{20},Din{21},Din{22},Din{23},Din{24},Din{25},Din{26},Din{27},Din{28},Din{29},Din{30},Din{31},Din{32},Din{33},Din{34},Din{35},Din{36},Din{37},Din{38});
% to build genedata matrix
[bool1,edge1(:,1)] = ismember(Pro,Din{2});
Matrix=zeros(number,36);
for i=1:number
    if bool1(i,1)==1
       for j =1:36
            Matrix(i,j)=express(edge1(i,1),j);  % genedata matrix
       end
    end
end
standarddeviation=zeros(number,1);
standarddeviation=std(Matrix,0,2);   
meanvalue=zeros(number,1);
meanvalue=mean(Matrix,2);
[bool,edge(:,1)] = ismember(Df_r{1},Pro); 
[bool,edge(:,2)] = ismember(Df_r{2},Pro);
AdjMatrix = zeros(number); 
for i = 1:length(Df_r{1}) 
    AdjMatrix(edge(i,1),edge(i,2)) = 1;
    AdjMatrix(edge(i,2),edge(i,1)) = 1;
end
SparseAdjMatrix = sparse(AdjMatrix); 
PCCMatrix=zeros(number);
ppc=zeros(length(Df_r{1}),1);
for i = 1:length(Df_r{1}) 
    for j=1:36
        if standarddeviation(edge(i,1),1) ~=0 & standarddeviation(edge(i,2),1) ~=0
            ppc(i,1) =  ppc(i,1)+(( Matrix(edge(i,1),j)-meanvalue(edge(i,1),1))/standarddeviation(edge(i,1),1) * ( Matrix(edge(i,2),j)-meanvalue(edge(i,2),1))/standarddeviation(edge(i,2),1)) ;
        end
    end 
    PCCMatrix(edge(i,1),edge(i,2)) = ppc(i,1)/35; 
    PCCMatrix(edge(i,2),edge(i,1)) = ppc(i,1)/35;
end
sumpcc=sum(PCCMatrix,2);
degree = sum(AdjMatrix,2); 
ECC= zeros(number,number); 
for i = 1:number
    for j = 1:number
         D=[Matrix(i,:);Matrix(j,:)];
        Jaccard(i,j)=1-pdist(D,'Jaccard');
    neighbornumber=0;
     if ( AdjMatrix(i,j) ~=0)
      for k=1:number
       if(k~=j)&&(k~=i)&&(AdjMatrix(i,k)~=0)&&(AdjMatrix(j,k)~=0)
       neighbornumber=neighbornumber+1;
       end
      end
      if((degree(i)>1) && (degree(j)>1))
      ECC(i,j)=neighbornumber/(sum(AdjMatrix(i,:))+sum(AdjMatrix(j,:))-2-neighbornumber);
      end
     end 
    end
end
SOECC = sum(ECC,2); 
Jaccard(isnan(Jaccard)==1)=0;
Jaccard=sum(Jaccard,2);
orthology= zeros(number,1); 
lengthprotein=length(protein);
for j=1:lengthprotein
    for i=1:number   
        if strcmp(protein{j,2},Pro{i})
            
            orthology(i,1)=str2num(protein{j,3})/99;
            continue;
        end
    end    
end
ESpro = importdata('Essential_yeast.txt');
filter =[100,200,300,400,500,600]; 
Jaccard=Jaccard/max(Jaccard);
sumpcc=sumpcc/max(sumpcc);
SOECC=SOECC/max(SOECC);
tp1= zeros(number,1);
for i=1:number
    if sumpcc(i)>SOECC(i)
        a=sumpcc(i);
    else
        a=SOECC(i);
    end
    
    if a>Jaccard(i)
       tp1(i) =Jaccard(i);
    else
        tp1(i)=a;
    end
end
s=0.05;
RWPOCPF=zeros(number,11);
for k=0:10
   for i=1:number
      for j=i+1:number
           if tp1(i,1)>tp1(j,1)-s && orthology(i,1)>orthology(j,1)-s 
              RWPOCPF(i,k+1)=RWPOCPF(i,k+1)+(tp1(i,1)-tp1(j,1))+(orthology(i,1)-orthology(j,1));
              
           end
           if tp1(i,1)-s<tp1(j,1) && orthology(i,1)-s<orthology(j,1) 
              RWPOCPF(j,k+1)=RWPOCPF(j,k+1)+(tp1(j,1)-tp1(i,1))+(orthology(j,1)-orthology(i,1));
           end
      end     
  end
end
value1 = zeros(number,11); 
inx1 = zeros(number,11);
RWPOCESnum = zeros(6,11);

for i=1:11
  [value1(:,i),inx1(:,i)] = sort(RWPOCPF(:,i),'descend'); %按分值按降序排列
end

for j=1:11
  for i = 1:6
     
     RWPOCESnum(i,j) = sum(ismember(Pro(inx1(1:filter(i),j)),ESpro)); 
     
  end
end



