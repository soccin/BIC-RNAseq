# library(testthat)

context("cut_lower_fun works")

# library(dendextend)
# library(dendextendRcpp)

test_that("cut_lower_fun works",{
   
      dend <- datasets::iris[1:4,-5] %>% dist %>% hclust %>% as.dendrogram
      # this is really cool!
      expect_identical(
         cut_lower_fun(dend, .4, labels),
         lapply(cut(dend, h = .4)$lower, labels)   
      )
   
      expect_identical(
         cut_lower_fun(dend, 0, labels),
         lapply(cut(dend, h = 0)$lower, labels)   
      )
   
      if("package:dendextendRcpp" %in% search()) {
      
   	   # cut should have returned the tree itself
   	   # but it FORCES a cut - as opposed to the cut_lower_fun (within the dendextendRcpp) function...
   	   expect_false(
   		  identical(
   			 dendextendRcpp::dendextendRcpp_cut_lower_fun(dend, 40, labels),
   			 lapply(cut(dend, h = 40)$lower, labels)   
   		  )
   	   )
   	   
   	   # cut should have returned the tree itself - instead it returns list()
   	   # as opposed to the cut_lower_fun (in the dendextendRcpp package) function...
   	   expect_false(
   		  identical(
   			 dendextendRcpp::dendextendRcpp_cut_lower_fun(dend, -1, labels),
   			 lapply(cut(dend, h = -1)$lower, labels)   )
   	   )
   	   
         tmp <- dendextendRcpp::dendextendRcpp_cut_lower_fun(dend[[1]], .4, function(x) {x})
   	   expect_identical(tmp,
   		  list(dend[[1]])  
   	   )

         # Here we should get:
            #  In Rcpp_get_dend_heights(tree, branches_heights = TRUE, labels_heights = FALSE) :
            #    	'height' attr is missing from node, 0 is returned, please check your tree.
#          expect_warning(all.equal(tmp,
#                    list(dend[[1]])))
         
   	}  # run only if dendextendRcpp was loaded

      
      # returns itself as it should:
      expect_identical(
         dendextend::dendextend_cut_lower_fun(dend[[1]], .4, function(x)x),
         list(dend[[1]])  
      )
      # this is the way to test this:
   #    expect_identical(
   #       old_cut_lower_fun(dend[[1]], .4, function(x)x),
   #       list(dend[[1]])  
   #    )   
      
      # function on a leaf gives what we expect
      expect_identical(
         cut_lower_fun(dend[[1]], .4, labels),
         list("1"))
      
      # cut, however, returns an empty list instead of itself...
      expect_identical(
         cut(dend[[1]], h = .4)$lower,
         list())
   
      
      
      
      #    library(microbenchmark)
   #    microbenchmark(
   #       cut_lower_fun(dend, .4, order.dendrogram),
   #       lapply(cut(dend, h = .4)$lower, order.dendrogram)   
   #    )
      #    # Rcpp is 3 times faster...
   #    
   
      
   #    dend_big = as.dendrogram(hclust(dist(iris[1:150,-5])))
   #       microbenchmark(
   #          cut_lower_fun(dend_big, .4, order.dendrogram),
   #          lapply(cut(dend_big, h = .4)$lower, order.dendrogram)   
   #       )
         # Rcpp is 16 times faster...
   
   
   
})





test_that("cut_lower_fun in dendextend",{
   
   dend <- as.dendrogram(hclust(dist(c(1,1,1,2,2))))
   expect_identical(dendextend_cut_lower_fun(dend, -.5, labels),
                    list())
   # But this will NOT be the same with dendextendRcpp !!
   
})