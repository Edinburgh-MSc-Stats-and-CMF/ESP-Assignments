setwd("/Users/share_xy/Desktop/Learning Materials/ESP-Assignments")
a <- scan("4300-0.txt", what="character", skip=73, nlines=32858-73, 
          fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("
split_punct <- function(text){
  p <- a[grep("[,.;!:?]", a)] ## return all words containing a punctuation mark
  words <- gsub("[,.;!:?]", "", p) ## get words without punctuation
  p <- gsub("[^,.;!:?]", "", p) ##retain the punctuation in the original word
  pwords <- paste(words, p) ## forming word and punctuation combinations
  a[grep("[,.;!:?]", a)] <- pwords ## reinsert the word and punctuation into a
  a <- strsplit(a, " ") ## split words by space and return a list
  a <- unlist(a, use.names = FALSE) ## Flatten the list to a vector
  a <- a[!(a %in% "")] #remove empty strings
}

a <- split_punct(a)
a

lowa <- tolower(a)
uniq_a <- unique(lowa)
uniq_a
m_uniq <- match(lowa, uniq_a)
m_uniq
tabul_uniq <- tabulate(m_uniq)

count <- table(lowa)
count <- count[order(count,decreasing = TRUE)]
b <- names(count)[1:1000]
#参数设置
mlag= 5
generated_num = 100

#生成matrix
M <- matrix(NA, nrow = length(a) - mlag, ncol = mlag + 1)
index <- match(lowa,b)
index
M <- t(sapply(1:(length(index) - mlag), function(i) index[i:(i + mlag)]))
## 把前面这个向量依次切片，生成一个向量，再转置
## M矩阵是单词的索引值（1-1000或NA)构成的矩阵

#初始化
generated_text = NULL
initial_row <- sample(which(!is.na(M[,1])), 1) ## 随机选择第一列不是NA的索引
generated_text[1] <- b[M[initial_row, 1]]  #随机选择一个索引1-1000的单词作为初始单词cd

for (i in 2:generated_num ) {
  # 标记是否生成了下一个单词
  word_generated <- FALSE
  
  # 从最高阶逐步退回到低阶
  for (j in mlag:1) {
    if (i > j) {  # 跳过滞后过长的情况
      # 提取前 j 个单词的索引
      prev_word_indices <- match(generated_text[(i-j):(i-1)], b)
      
      # 在矩阵 M 中查找匹配的行
      matching_rows <- which(apply(M[, 1:j, drop = FALSE], 1, function(x) all(x == prev_word_indices)))
      
      # 如果找到匹配的行，从中随机选择下一个单词
      if (length(matching_rows) > 0) {
        next_word_index <- sample(M[matching_rows, j + 1], 1)
        
        # 检查是否有非 NA 值
        if (!is.na(next_word_index)) {
          generated_text[i] <- b[next_word_index]
          word_generated <- TRUE
          break  # 成功生成单词，退出 j 循环
        }
      }
    }
  }
  # 如果所有阶都没有成功生成单词，随机选择一个非 NA 单词
  if (!word_generated) {
    non_na_indices <- which(!is.na(M[, 1]))
    random_row <- sample(non_na_indices, 1)
    generated_text[i] <- b[M[random_row, 1]]
  }
}

sentence <- paste(generated_text, collapse = " ")
# 去除标点符号前的空格，例如 " , " 替换为 ","
sentence <- gsub(" ([,.!?;:])", "\\1", sentence)
print(sentence)