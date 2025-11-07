#help function1:calculate length for S ended by end_num
calculate_length <- function(end_num) {
  if (end_num <= 0) return(0)
  k <- floor(log10(end_num)) + 1
  
  # length of S:S is constructed by all k digit number
  if (k > 1) {
    i <- 1:(k-1)
    complete_part <- sum(9 * (10^(i-1)) * i)
  } else {
    complete_part <- 0
  }
  
  incomplete_part <- (end_num - 10^(k-1) + 1) * k
  
  return(complete_part + incomplete_part)
}

#help function2:check if first appears between 1-9  
check_1digit <- function(s,n){
  # Case 1: Single digit (1-9)
  if (n == 1) {
    digit <- as.integer(s)
    if (digit >= 1 && digit <= 9) {
      return(digit)
    }
    if (digit==0){
      return(11)
    }
  }
  if (n>1){
      first_num_str <- substr(s, 1, 1)
      #print(first_num_str)
      first_num <- as.integer(first_num_str)
      concat <- ""
      num <- first_num
      # 拼接连续数字，直到拼接长度 >= n
      while (nchar(concat) < n) {
        concat <- paste0(concat, num)
        num <- num + 1
        #print(concat)
      }
      if (substr(concat, 1, n) == s) {
        return(first_num)
      }
    }
    return(FALSE)
  }

#help function3:find sort of compare
find_test_sort <- function(s) {
  digits <- as.integer(strsplit(s, "")[[1]])
  positions <- seq_along(digits)
  nonzero_index <- which(digits != 0)
  digits_nonzero <- digits[nonzero_index]
  pos_nonzero <- positions[nonzero_index]
  sorted_index <- order(digits_nonzero, pos_nonzero)
  sorted_positions <- pos_nonzero[sorted_index]
  return(sorted_positions)
}

#help function4:向后补位
add_number <- function(s,n,len_check,start_pos){
  sub_s <- substr(s, start_pos, n)
  len_diff <- len_check - nchar(sub_s)
  if (len_diff > 0) {
    # 获取要补的数字
    prev_num_str <- substr(s,start_pos - len_diff,start_pos - 1)
    add_num <- as.character(as.integer(prev_num_str) + 1)
    # 如果位数比差的大，截取后diff位（遇9取0）
    if (nchar(add_num) > len_diff) {
      pad_num <- substr(add_num,nchar(add_num)-len_diff+1, nchar(add_num))
    }else if (nchar(add_num) < len_diff){
      pad_num <- sprintf("%0*d", len_diff, as.integer(add_num))
    }else {
      pad_num <- add_num
    }
      sub_s <- paste0(sub_s,pad_num)
  }
  num_full <- as.integer(sub_s)
  return(num_full)
}

# help function6:向前比较

success_forward <- function(start_number,s,check_posi,len_check){
  offset <- 0
  success_forward <- TRUE
  curr_chunk <- start_number
  #print(curr_chunk)
  
  if (log10(as.integer(curr_chunk)) %% 1 == 0){
    check_start <- check_posi - offset - len_check+1
    check_end <- check_posi - offset - 1
    if (check_end < 1) {
      return(success_forward)
    }
    prev_chunk <- substr(s, max(1, check_start), check_end)
    #print(prev_chunk)
    if (nchar(prev_chunk) <= 0) {
      return(success_forward)
    }
    
    forward_check_length <- nchar(prev_chunk)
    last_digit <- as.integer(substr(curr_chunk,len_check - forward_check_length + 1,len_check))
    prev_digit <- as.integer(prev_chunk)
    #print(last_digit)
    if(last_digit==0){
      last_digit0 <- substr(curr_chunk,len_check - forward_check_length + 1,len_check)
      last_digit <- as.integer(paste0("1", last_digit0))
      #print(last_digit)
    }
    if (!(prev_digit == (last_digit - 1) || (last_digit == 0 && prev_digit == 9))) {
      success_forward <- FALSE
      return(success_forward)
    }
    len_check <- len_check-1
    offset <- offset + len_check
    curr_chunk <- prev_chunk
  }
  
  # 只要没到头就继续比
  while (TRUE) {
  
    check_start <- check_posi - offset - len_check
    #print(len_check)
    check_end <- check_posi - offset - 1
    #print(check_start)
    #print(check_end)
    if (check_end < 1) break
    
    prev_chunk <- substr(s, max(1, check_start), check_end)
    #print(prev_chunk)
    # print(curr_chunk)
    if (nchar(prev_chunk) <= 0) break
    
    forward_check_length <- nchar(prev_chunk)
    last_digit <- as.integer(substr(curr_chunk,len_check - forward_check_length + 1,len_check))
    prev_digit <- as.integer(prev_chunk)
    # print(last_digit)
    # print(last_digit)
    if(last_digit==0){
      last_digit0 <- substr(curr_chunk,len_check - forward_check_length + 1,len_check)
      last_digit <- as.integer(paste0("1", last_digit0))
      #print(last_digit)
    }
    
    # 检查递减规律
    if (!(prev_digit == (last_digit - 1) ||(last_digit == 0 && prev_digit == 9))) {
      success_forward <- FALSE
      break
    }
    offset <- offset + len_check
    curr_chunk <- prev_chunk
  }
  #print(len_check)
  return(success_forward)
}


# help function7:向后比较
success_backward <- function(start_number, s, check_posi, len_check) {
  offset <- 0
  success_backward <- TRUE
  curr_chunk <- start_number
  n <- nchar(s)
  
  while (TRUE) {
    next_start <- check_posi + offset + len_check
    next_end <- next_start + len_check - 1
    if (next_start > n) break
    
    next_chunk <- substr(s, next_start, min(next_end, n))

    if (nchar(next_chunk) < len_check) {
      num_full <- add_number(s, n, len_check, next_start)
      if (num_full == (as.integer(curr_chunk) + 1)) {
        success_backward <- TRUE
      } else {
        success_backward <- FALSE
      }
      break
    }
    if (as.integer(next_chunk) != (as.integer(curr_chunk) + 1)) {
      success_backward <- FALSE
      break
    }

    offset <- offset + len_check
    curr_chunk <- next_chunk
  }
  
  return(success_backward)
}



num_str <- function(s) {
  n <- nchar(s)
  if (n == 0) return(NA)
  if (n>0){
    check_status <- check_1digit(s,n)
    if (check_status!=FALSE){
      return(check_status)
    }
    
    if (check_status==FALSE){
      for (len_check in 2:n){
        # check 从头开始是否有规律
        first_number <- as.integer(substr(s, 1, 1+len_check-1))
        
        # 如果check的长度和s的长度一样
        if (len_check==n){
          test_sort <- find_test_sort(s)
          if (length(test_sort) == 0){
            num <- as.integer(paste0("1", s))
            num_posi <- calculate_length(num)
            final_posi_special <- num_posi-len_check+1
            return(final_posi_special)
          }

          if (as.integer(substr(s,1,1))!=0){
          num_posi <- calculate_length(as.integer(s))
          final_posi_first_posi <- num_posi-len_check+1
          }
          
          if (1 %in% test_sort) {
            test_sort <- test_sort[test_sort != 1]
          }
          
          if (length(test_sort) == 0){
            return(final_posi_first_posi)
          }else{
            for (check_posi in test_sort){
              first_chunk <- substr(s, check_posi, check_posi + len_check - 1)
              first_num <- as.integer(first_chunk)
              
              #从第一个chunk开始就补位
              if (nchar(first_chunk) < len_check){
                num_full <- add_number(s,n,len_check,start_pos=check_posi)
                #print(num_full)
                #向前检查
                forword_status <- success_forward(start_number =  as.character(num_full),s,check_posi=check_posi,len_check=len_check)
                #print(forword_status)
                if (forword_status){
                  num_posi <- calculate_length(num_full)
                  final_posi <- num_posi - (check_posi-1) - len_check + 1
                  if (as.integer(substr(s,1,1))!=0){
                    return(min(final_posi_first_posi,final_posi))
                  }else{
                    return(final_posi)
                  }
                }else{
                  next  # 进入下一次循环
                }
              }
            }
          }
        }
        
        # 如果check的长度小于s的长度
        if (len_check < n){
          test_sort <- find_test_sort(s)
          for (check_posi in test_sort){
            first_chunk <- substr(s, check_posi, check_posi + len_check - 1)
            first_num <- as.integer(first_chunk)
            
            #从第一个chunk开始就补位
            if (nchar(first_chunk) < len_check){
              num_full <- add_number(s,n,len_check,start_pos=check_posi)
              #print(num_full)
              #print(len_check)
              #向前检查
              forword_status <- success_forward(start_number =  as.character(num_full),s,check_posi=check_posi,len_check=len_check)
              #print(forword_status)
              if (forword_status){
                num_posi <- calculate_length(num_full)
                final_posi <- num_posi - (check_posi-1) - len_check + 1
                return(final_posi)
              }else{
                next  # 进入下一次循环
              }
            }
            
            if (nchar(first_chunk) == len_check){
              #向后检查
              backward_status <- success_backward(start_number =  first_chunk, s,check_posi,len_check)
              forword_status <- success_forward(start_number =  first_chunk,s,check_posi=check_posi,len_check=len_check)
              if(backward_status&forword_status){
                num_posi <- calculate_length(as.integer(first_chunk))
                final_posi <- num_posi - (check_posi-1) - len_check + 1
                return(final_posi)
              }else{
                next
              }
            }
          
          }
        }
     
      }
    }
  }
}

# test_cases <- c(
#   "1","5","12","1234","78",
#   "910","1011","99100","9899100",
#   "517617","645","702","299300",
#   "4321","210","321",
#   "99","100","100101","9991000","1000999",
#   "12345678910","9101112","987654","909192","89012","449","99200",
#   "1", "12", "123", "1234", "12345",
#   # 含零
#   "102", "103", "110", "112", "130",
#   "100", "101", "1000", "1001", "1011",
#   # 多个1
#   "111", "1111", "1112",
#   "001", "0", "00", "002", "010",
#   # 含零
#   "46513", "431", "91011", "8991", "20321"
# )
# 
# for (s in test_cases) {
#   cat("num_str(\"", s, "\") = ", num_str(s), "\n", sep = "")
# }


