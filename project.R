
rm(list=ls())
library(utils)
library(tibble)
library(magrittr)
library(tidyverse)
library(tm)
library(wordcloud)
library(ggplot2)
library("scales")
#reading the csv files
nmax=100000000
train_variants<- read.csv("training_variants.csv")
training_text <- read_delim("training_text.csv", delim = '|', 
                            skip = 1, col_names = c("ID", "VV", "Text"),
                            n_max = nmax)
training_text$VV = NULL
test_variants <- read.csv("test_variants.csv")
test_text <- read_delim("test_text.csv", delim = '|', 
                        skip = 1, col_names = c("ID", "VV", "Text"), 
                        n_max = nmax)
test_text$VV = NULL
#Exploratory Data analysis
str(train_variants)
colnames(train_variants)
head(train_variants)

str(training_text)
colnames(training_text)
head(training_text)

str(test_variants)
colnames(test_variants)
head(test_variants)

str(test_text)
colnames(test_text)
head(test_text)

ggplot(train_variants,aes(x=Class)) + geom_bar()

ggplot(train_variants, aes_string(x = train_variants$Gene, y = train_variants$Variation)) + 
  geom_point(aes_string(colour = train_variants$Class),size=4) +
  theme_bw()+ ylab("Variants") + xlab("Gene") + ggtitle("Scatter") + 
  theme(text=element_text(size=25)) + 
  scale_x_continuous(breaks=pretty_breaks(n=10)) +
  scale_y_continuous(breaks=pretty_breaks(n=10)) +
  scale_colour_discrete(name="Class") + 
  scale_shape_discrete(name="CAD")

#feature engineering 

train_variants$ID<-as.factor(train_variants$ID)
train_variants$Class<-as.factor(train_variants$Class)
train_variants$Gene<-as.numeric(train_variants$Gene)
train_variants$Variation<-as.numeric(train_variants$Variation)

test_variants$ID<-as.factor(test_variants$ID)
test_variants$Variation<-as.numeric(test_variants$Variation)
test_variants$Gene<-as.numeric(test_variants$Gene)
test_variants$Class = NA
test_variants$Class = as.factor(test_variants$Class)

#combining training text and test text
txt <- rbind(training_text, test_text)

#storing row count of text file
traincnt<-nrow(training_text)
testcnt <- nrow(test_text)

#text mining
tcorpus = Corpus(VectorSource(txt$Text))
tcorpus = tm_map(tcorpus, stripWhitespace)
tcorpus = tm_map(tcorpus, removePunctuation)
tcorpus = tm_map(tcorpus, removeNumbers)
tcorpus = tm_map(tcorpus, content_transformer(tolower))
tcorpus = tm_map(tcorpus, removeWords, stopwords('english'))
tcorpus = tm_map(tcorpus, stemDocument, language="english")

#wordcloud
graphics.off()
wordcloud(tcorpus, max.words = 200, scale=c(3, .1), colors=brewer.pal(6, "Dark2"))



#Converting the corpus into document term matrix
dtm <- DocumentTermMatrix(tcorpus, control = list(weighting = weightTfIdf))
dtms <- removeSparseTerms(dtm, 0.99)
dat_dtms <- as.matrix(dtms)

#spliting output into train and test
trainout <- dat_dtms[1:traincnt, ]
testout <- dat_dtms[(traincnt + 1):(traincnt + testcnt), ]

#combining Documenttermmatrix with variants file
train <- cbind(train_variants, trainout)
test <- cbind(test_variants, testout)


train$ID = NULL
test$ID = NULL
#storing class list
cl=unique(train$Class)

#dividing train to trainmodel and testmodel
library(caTools)
set.seed(101) 
sample = sample.split(train$Class, SplitRatio = 0.7)
trainmodel = subset(train, sample == TRUE)
testmodel  = subset(train, sample == FALSE)

library(e1071)
#naive bayes 
model = naiveBayes(Class ~ ., data = trainmodel,type="prob")
#predict using test data
test_pred = predict(model, testmodel[,-3], type = "class")
test_pred_prob = predict(model, testmodel[,-3], type = "prob")


#Visualizing the confusion matrix
xtab = table(observed = testmodel[,3], predicted = test_pred)#u can predict accuracy only on class not on probabilities
library(caret)
confusionMatrix(xtab)

#original test data prediction
predText_prob = predict(model, test[,-3], type = "raw")
#saving output into file
pred <- matrix(predText_prob , ncol = length(cl), byrow = TRUE)
subfil<- cbind(test_variants$ID, pred)
str(subfil)
colnames(subfil)<-c("ID"	,"class1",	"class2",	"class3",	"class4","class5","class6",	"class7",	"class8","class9")
write.csv(subfil,"Outputfile.csv",row.names = F)

