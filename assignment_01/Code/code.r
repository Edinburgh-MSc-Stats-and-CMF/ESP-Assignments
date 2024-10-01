setwd("/Users/share_xy/Desktop/Learning Materials/ESP-Assignments")
a <- scan("4300-0.txt", what="character", skip=73, nlines=32858-73, 
          fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) # remove "_("

## write a function to separate punctuation from words in text
split_punct <- function(text){
  p <- a[grep("[,.;!:?]", a)] # return words with punctuation 
  words <- gsub("[,.;!:?]", "", p) # remove punctuation from words
  p <- gsub("[^,.;!:?]", "", p) # retain only punctuation
  pwords <- paste(words, p) # recombine words with punctuation(with spaces)
  a[grep("[,.;!:?]", a)] <- pwords # replace original ones with word-punctuation forms
  a <- strsplit(a, " ") # split words by spaces to separate punctuation
  a <- unlist(a) # flatten into a vector
  a <- a[!(a %in% "")] # remove empty strings
}
a <- split_punct(a) # apply split_punct to text

## select the most frenquent words from text
lowa <- tolower(a) # convert text to lowercase
uniq_a <- unique(lowa) # create a unique word list
m_uniq <- match(lowa, uniq_a) # map each word in 'lowa' to its index in 'uniq_a'
tabul_uniq <- tabulate(m_uniq) # tabulate the frequency of each unique word

count <- table(lowa) # create a frequency table
count <- count[order(count,decreasing = TRUE)] # sort it in decreasing order
b <- names(count)[1:1000] # select the top 1000 most frenquent words

## set parameters for Markov chain text generation
mlag <- 6  # number of previous words to consider(Markov lag)
generated_num <- 100 # number of words to generate

## make the matrices of common word token sequences
M <- matrix(NA, nrow = length(a) - mlag, ncol = mlag + 1) # initialize an empty matrix M
index <- match(lowa,b) # map each word in 'lowa' to its index in 'b'
M <- t(sapply(1:(length(index) - mlag), function(i) index[i:(i + mlag)])) # create a transposed matrix M with overlapping sequences of 'index' of length mlag + 1
M <- M[!is.na(M[,1]),] # remove rows with NA in column 1


## generate random text
generated_text <- "he" # start with the word "he"

for (i in 2:generated_num ) {
  word_generated <- FALSE # tract whether a word has been generated
  # iterate from higher-order Markov models to lower-order ones
  for (j in mlag:1) { 
    # skip lags that are too long
    if (i > j) {  
      prev_word_indices <- match(generated_text[(i-j):(i-1)], b) # extract the index of the first j words
      # find rows in M that match the previous word sequence
      matching_rows <- M
      for (k in 1:j) {
        matching_rows <- matching_rows[matching_rows[,k] %in% prev_word_indices[k],] # narrow down matching rows
        matching_rows <- matrix(matching_rows,ncol = ncol(M)) # retain matrix structure
      }
      matching_rows <- matching_rows[!is.na(matching_rows[,j+1]),]  # retain rows where the next word is non-NA
      matching_rows <- matrix(matching_rows,ncol = ncol(M)) # retain matrix structure
      
      # if matching rows are found, randomly select the next word from them
      if (length(matching_rows) > 0) {
        next_word_index <- sample(matching_rows[, j + 1], 1)
        # check for non-NA values
        if (!is.na(next_word_index)) {
          generated_text[i] <- b[next_word_index]
          word_generated <- TRUE
          print(j) # debugging output to track the order of Markov chain
          break  # exit the loop once a word is generated
        }
      }
    }
  }
  # if no matching sequence is found, randomly select a word from M
  if (!word_generated) {
    random_row <- sample(M[,1], 1)
    print("rm") # debugging output indicating random word selection
    generated_text[i] <- b[M[random_row, 1]]
  }
}
sentence <- paste(generated_text, collapse = " ") # combine the generated words into a sentence
sentence <- gsub(" ([,.!?;:])", "\\1", sentence) # removing extra spaces before punctuation
print(sentence) # print the generated sentence

