slots <- slotNames(x)
for (slot_name in slots) {
  slot_content <- slot(x, slot_name)
  cat(paste0("--- Slot: ", slot_name, " ---
"))
  cat(paste0("Class: ", class(slot_content), "\n"))
  if (is(slot_content, "matrix") || is(slot_content, "array") || is(slot_content, "data.frame")) {
    cat(paste0("Dimensions: ", paste(dim(slot_content), collapse = " x "), "\n"))
    # Preview the head of tabular data if available
    if ((is(slot_content, "data.frame") || is(slot_content, "matrix")) && nrow(slot_content) > 0 && ncol(slot_content) > 0) {
      cat("Head:\n")
      print(head(slot_content))
    } else if (nrow(slot_content) == 0 || ncol(slot_content) == 0) {
      cat("Object is empty.\n")
    }
  } else if (is(slot_content, "list")) {
    cat(paste0("Number of elements: ", length(slot_content), "\n"))
  } else {
    cat("Content summary:\n")
    print(slot_content)
  }
  cat("\n") # Add a newline for better separation
}
