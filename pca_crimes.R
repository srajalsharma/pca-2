library(factoextra)
p=21
a=c(122679, 129017, 583645, 938110, 980858, 897171, 989071, 1040046, 1050945)
b=c(77564, 86469, 81483, 76248, 76555, 72829, 78051, 76851, 63359)
c1=c(228650, 244270, 309546, 339457, 329243, 338954, 359849, 378236, 405861)
d=c(33098, 38172, 58224, 89423, 94172, 106958, 129032, 141764, 148185)
e=c(25125, 27936, 31725, 38455, 33433, 35849, 33606, 31591, 32235)
f=c(22571, 22571, 22571, 18714, 20532, 21410, 22727, 24349, 27696)
g=c(33719, 39577, 15365, 40401, 38670, 40801, 1400, 1340, 1670)
h=c(148083, 148083, 129306, 142560, 150170, 143524, 148972, 156268, 165782)
i=c(4514, 4514, 4514, 4966, 5250, 4439, 4062, 4129, 4243)
j=c(2213, 3477, 5693, 9622, 11592, 12317, 21796, 27248, 44546)
k=c(7256, 7256, 7256, 5396, 6040, 6986, 9013, 8536, 7569)
l=c(21288, 21288, 21288, 5835, 5156, 4732, 42143, 35196, 34671)
m=c(1766, 1766, 1766, 1473, 1278, 1226, 2242, 2128, 2251)
n1=c(3517, 3554, 3940, 5466, 6877, 8132, 2854, 2278, 2269)
o=c(378843, 378843, 378843, 78843, 78843, 660577, 368616, 414658, 71523)
p1=c(31880, 34756, 516648, 828495, 828495, 796032, 775263, 802372, 854618)
q=c(106141, 106141, 106141, 106141, 106141, 15205, 128601, 135812, 144949)
r=c(71365, 71365, 71365, 71365, 71365, 53929, 92107, 66305, 73122)
s=c(424129, 424129, 424129, 59098, 57668, 60560, 424129, 1266437, 676883)
t=c(295617, 295617, 295617, 295617, 295617, 281830, 355994, 257243, 287404)
u=c(270548, 270548, 270548, 700046, 757622, 25147, 3158, 450, 867)

#scaling of data
data=matrix(c(a,b,c1,d,e,f,g,h,i,j,k,l,m,n1,o,p1,q,r,s,t,u), nrow = 9, ncol = 21, byrow = FALSE)
y=abs(scale(data))

#covariance matrix
co=cov(y)

#testing: lawley procedure
cm=cor(y)
i=1
ri=c(1:21)*0
for(i in 1:21)
{
ri[i]=c((sum(cm[,i])-1)/(p-1)) #ri
}
rj=(sum(cm)-21)/(p*(p-1))
rj #r_bar
s2=0
for(i in 1:21)
{
for(j in 1:21)
{
if(i!=j)
{
s2=s2+((cm[i,j]-rj)*(cm[i,j]-rj))
}
}
}
s3=s2/2   #summation 1
s4=0
for(i in 1:21)
{
s4=s4+((ri[i]-rj)^2)   #summation 2
}
gaa=(((p-1)^2)*(1-(1-rj)^2))/(p-(p-2)*((1-rj)^2))  #gamma
n=9
T=((n-1)/((1-rj)^2))*(s3-(gaa*s4))
T
T_tab=243.727  #chi-square value at 5% los with ((p+1)(p-2))/2=209 dof 

if(T>T_tab)
{
print(' we reject our null hypothesis, hence all the eigen values are different and pca is applicable')
} else { print('pca cannot be applied')
}

#eigenvalues and eigenvectors
ei=eigen(co)
eva=ei$values
eve=ei$vectors

#principal components
result=t(eve)%*%t(y)

#variation explained by principal components
pc1.pca=prcomp(y, scale=FALSE)
summary(pc1.pca)
fviz_eig(pc1.pca)

#principal components to be chosen
p=21
print('Principal Components that effectively summarises total variance are:')
for(v in 1:p)
{
if(round(eva[v],1)>round(mean(eva),1))
{
print(v)
w=v
}
}
print('Their value is represented by each column of following matrix:')
result[,1:w]

#correlation between original variables x and principal components y
cpc=matrix(nrow=21,ncol=9,byrow=TRUE,dimnames=list(c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","y16","y17","y18","y19","y20","y21"),c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9")))
for(i in 1:21)
{
for(j in 1:9)
{
cpc[i,j]=((sqrt(eva[j]))*(eve[i,j]))/(sqrt(co[i,i]))
}
}
print(cpc)

