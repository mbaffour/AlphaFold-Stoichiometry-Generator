# ══════════════════════════════════════════════════════════════════════════════
# AlphaFold Stoichiometry Generator
# Build combinatorial AlphaFold multimer inputs with full control
# ══════════════════════════════════════════════════════════════════════════════

library(shiny)
library(DT)
library(jsonlite)
library(zip)

# ── Helper Functions ──────────────────────────────────────────────────────────

parse_fasta <- function(text) {
  lines <- strsplit(trimws(text), "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0]
  if (length(lines) == 0) return(NULL)
  
  proteins <- list()
  current_id <- NULL
  current_seq <- ""
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.null(current_id)) {
        proteins[[length(proteins) + 1]] <- list(id = current_id, sequence = current_seq)
      }
      current_id <- sub("^>\\s*", "", line)
      current_id <- sub("\\s.*", "", current_id)
      current_seq <- ""
    } else {
      clean <- toupper(gsub("[^A-Za-z]", "", line))
      current_seq <- paste0(current_seq, clean)
    }
  }
  
  if (is.null(current_id) && nchar(paste0(lines, collapse = "")) > 0) {
    seq_text <- toupper(gsub("[^A-Za-z]", "", paste0(lines, collapse = "")))
    proteins[[1]] <- list(id = "Protein_1", sequence = seq_text)
  } else if (!is.null(current_id)) {
    proteins[[length(proteins) + 1]] <- list(id = current_id, sequence = current_seq)
  }
  
  if (length(proteins) == 0) return(NULL)
  
  # Deduplicate IDs
  ids <- sapply(proteins, function(p) p$id)
  if (any(duplicated(ids))) {
    for (i in seq_along(proteins)) {
      dupes <- which(ids == ids[i])
      if (length(dupes) > 1) {
        for (j in seq_along(dupes)) {
          proteins[[dupes[j]]]$id <- paste0(ids[dupes[j]], "_", j)
        }
        ids <- sapply(proteins, function(p) p$id)
      }
    }
  }
  
  proteins
}

validate_sequence <- function(seq) {
  valid_aa <- "ACDEFGHIKLMNPQRSTVWY"
  invalid <- gsub(paste0("[", valid_aa, "]"), "", toupper(seq))
  if (nchar(invalid) > 0) {
    unique_bad <- unique(strsplit(invalid, "")[[1]])
    return(paste("Invalid characters:", paste(unique_bad, collapse = ", ")))
  }
  return(NULL)
}

# Build readable ratio: "ProtA x 2 | ProtB x 1"
format_ratio_display <- function(counts, protein_ids) {
  parts <- c()
  for (i in seq_along(counts)) {
    if (counts[i] > 0) {
      parts <- c(parts, paste0(protein_ids[i], " x ", counts[i]))
    }
  }
  paste(parts, collapse = "  |  ")
}

# Build compact model name: "A2_B1"
format_model_name <- function(counts, protein_ids) {
  parts <- c()
  for (i in seq_along(counts)) {
    if (counts[i] > 0) {
      parts <- c(parts, paste0(protein_ids[i], counts[i]))
    }
  }
  paste(parts, collapse = "_")
}

generate_combinations <- function(proteins, ranges) {
  grid <- expand.grid(ranges)
  prot_ids <- sapply(proteins, function(p) p$id)
  colnames(grid) <- prot_ids
  
  combos <- list()
  for (i in seq_len(nrow(grid))) {
    row <- as.integer(grid[i, ])
    names(row) <- prot_ids
    total_res <- sum(row * sapply(proteins, function(p) nchar(p$sequence)))
    
    combos[[i]] <- list(
      model_name  = format_model_name(row, prot_ids),
      display     = format_ratio_display(row, prot_ids),
      counts      = row,
      total_chains = sum(row),
      total_residues = total_res,
      warning     = total_res > 3000
    )
  }
  combos
}

generate_random_combinations <- function(proteins, ranges, n = 10) {
  prot_ids <- sapply(proteins, function(p) p$id)
  combos <- list()
  seen <- c()
  attempts <- 0
  
  while (length(combos) < n && attempts < n * 20) {
    attempts <- attempts + 1
    row <- sapply(ranges, function(r) sample(r, 1))
    names(row) <- prot_ids
    sig <- paste(row, collapse = "-")
    if (sig %in% seen) next
    seen <- c(seen, sig)
    
    total_res <- sum(row * sapply(proteins, function(p) nchar(p$sequence)))
    combos[[length(combos) + 1]] <- list(
      model_name  = format_model_name(row, prot_ids),
      display     = format_ratio_display(row, prot_ids),
      counts      = row,
      total_chains = sum(row),
      total_residues = total_res,
      warning     = total_res > 3000
    )
  }
  combos
}

generate_custom_combinations <- function(proteins, ratio_text) {
  prot_ids <- sapply(proteins, function(p) p$id)
  n_prot <- length(proteins)
  lines <- strsplit(trimws(ratio_text), "\n")[[1]]
  lines <- trimws(lines[nchar(trimws(lines)) > 0])
  
  combos <- list()
  for (line in lines) {
    nums <- as.integer(strsplit(gsub("[^0-9]+", " ", trimws(line)), "\\s+")[[1]])
    nums <- nums[!is.na(nums)]
    if (length(nums) < n_prot) {
      nums <- c(nums, rep(0, n_prot - length(nums)))
    } else if (length(nums) > n_prot) {
      nums <- nums[seq_len(n_prot)]
    }
    names(nums) <- prot_ids
    total_res <- sum(nums * sapply(proteins, function(p) nchar(p$sequence)))
    
    combos[[length(combos) + 1]] <- list(
      model_name  = format_model_name(nums, prot_ids),
      display     = format_ratio_display(nums, prot_ids),
      counts      = nums,
      total_chains = sum(nums),
      total_residues = total_res,
      warning     = total_res > 3000
    )
  }
  combos
}

combo_to_fasta <- function(combo, proteins) {
  lines <- c()
  for (j in seq_along(proteins)) {
    prot <- proteins[[j]]
    count <- combo$counts[j]
    if (count == 0) next
    for (k in seq_len(count)) {
      lines <- c(lines, paste0(">", prot$id, "_chain_", k, " | Model: ", combo$model_name))
      seq_chars <- prot$sequence
      while (nchar(seq_chars) > 0) {
        lines <- c(lines, substr(seq_chars, 1, 80))
        seq_chars <- substr(seq_chars, 81, nchar(seq_chars))
      }
    }
  }
  paste(lines, collapse = "\n")
}

combo_to_json <- function(combo, proteins) {
  seq_grouped <- list()
  for (j in seq_along(proteins)) {
    prot <- proteins[[j]]
    count <- combo$counts[j]
    if (count > 0) {
      seq_grouped[[length(seq_grouped) + 1]] <- list(
        proteinChain = list(
          sequence = prot$sequence,
          count = count
        )
      )
    }
  }
  
  job <- list(
    name = combo$model_name,
    modelSeeds = list(1),
    sequences = seq_grouped
  )
  
  toJSON(list(list(job)), auto_unbox = TRUE, pretty = TRUE)
}

combo_to_summary <- function(combo, proteins) {
  lines <- c(
    paste("Model:", combo$model_name),
    paste("Stoichiometry:", combo$display),
    paste("Total chains:", combo$total_chains),
    paste("Total residues:", format(combo$total_residues, big.mark = ",")),
    ""
  )
  if (combo$warning) {
    lines <- c(lines, "** WARNING: Exceeds 3,000 residues -- AlphaFold may struggle **", "")
  }
  lines <- c(lines, "Chain Details:")
  lines <- c(lines, paste0(rep("-", 55), collapse = ""))
  for (j in seq_along(proteins)) {
    prot <- proteins[[j]]
    count <- combo$counts[j]
    if (count > 0) {
      lines <- c(lines, sprintf("  %-20s  %d copies  (%d aa each = %d aa total)",
                                prot$id, count, nchar(prot$sequence),
                                count * nchar(prot$sequence)))
    }
  }
  lines <- c(lines, paste0(rep("-", 55), collapse = ""))
  paste(lines, collapse = "\n")
}

# ── Custom CSS ────────────────────────────────────────────────────────────────

app_css <- "
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
  .app-header {
    background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
    color: white; padding: 18px 25px; margin: -15px -15px 20px -15px;
    border-radius: 0 0 8px 8px;
  }
  .app-header h2 { margin: 0 0 4px 0; font-weight: 700; font-size: 1.5em; letter-spacing: -0.5px; }
  .app-header p { margin: 0; opacity: 0.8; font-size: 0.9em; }
  .stat-card {
    background: #f8f9fa; border-left: 4px solid #0f3460; border-radius: 6px;
    padding: 12px 16px; margin-bottom: 10px;
  }
  .stat-card .stat-number { font-size: 1.8em; font-weight: 700; color: #0f3460; }
  .stat-card .stat-label { font-size: 0.8em; color: #666; text-transform: uppercase; letter-spacing: 0.5px; }
  .stat-card.warning { border-left-color: #e74c3c; }
  .stat-card.warning .stat-number { color: #e74c3c; }
  .stat-card.success { border-left-color: #27ae60; }
  .stat-card.success .stat-number { color: #27ae60; }
  .section-label {
    font-weight: 600; color: #0f3460; font-size: 1em;
    margin-bottom: 8px; padding-bottom: 5px; border-bottom: 2px solid #e8e8e8;
  }
  .sequence-badge {
    display: inline-block; background: #e8f4f8; color: #0f3460;
    padding: 3px 10px; border-radius: 12px; font-size: 0.85em;
    margin: 2px; font-weight: 500;
  }
  .btn-generate {
    background: linear-gradient(135deg, #27ae60, #2ecc71); border: none;
    color: white; font-weight: 600; letter-spacing: 0.3px;
  }
  .btn-generate:hover { background: linear-gradient(135deg, #219a52, #27ae60); color: white; }
  .btn-parse {
    background: linear-gradient(135deg, #0f3460, #16213e); border: none;
    color: white; font-weight: 600;
  }
  .btn-parse:hover { background: linear-gradient(135deg, #16213e, #1a1a2e); color: white; }
  .warning-text { color: #e74c3c; font-weight: 600; }
  .help-text { font-size: 0.82em; color: #888; margin-top: 3px; }
  .dl-section {
    background: #f0f7ff; padding: 16px; border-radius: 8px; margin-bottom: 15px;
    border: 1px solid #d0e3f7;
  }
  .protein-loaded {
    background: #d4edda; border: 1px solid #c3e6cb; padding: 10px 14px;
    border-radius: 6px; margin-bottom: 10px;
  }
  .stoich-protein-block {
    background: #f8f9fa; padding: 10px; border-radius: 6px;
    margin-bottom: 8px; border: 1px solid #e9ecef;
  }
  pre.shiny-text-output {
    background: #1e1e1e; color: #d4d4d4; border: none;
    border-radius: 8px; padding: 15px; font-size: 0.85em;
  }
"

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  tags$head(tags$style(HTML(app_css))),
  
  div(class = "app-header",
      h2("AlphaFold Stoichiometry Generator"),
      p("Design and export multimer modeling inputs with systematic stoichiometry control")
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # -- Input --
      div(class = "section-label", "Input Sequences"),
      radioButtons("input_mode", NULL,
                   choices = c("Paste FASTA" = "paste", "Upload File" = "upload"),
                   inline = TRUE),
      conditionalPanel(
        "input.input_mode == 'paste'",
        textAreaInput("fasta_text", NULL,
                      placeholder = ">GeneA\nMKFLILALF...\n>GeneB\nMSVQTIGDL...\n>GeneC\nMTDKLFPEN...",
                      rows = 10, width = "100%"),
        div(class = "help-text", "Paste FASTA or plain sequence. Multiple sequences supported.")
      ),
      conditionalPanel(
        "input.input_mode == 'upload'",
        fileInput("fasta_file", NULL, accept = c(".fasta", ".fa", ".faa", ".txt")),
        div(class = "help-text", "Upload .fasta, .fa, .faa, or .txt")
      ),
      actionButton("parse_btn", "Parse Sequences", class = "btn-parse btn-sm", width = "100%"),
      
      uiOutput("parsed_summary_sidebar"),
      
      hr(),
      
      # -- Stoichiometry --
      div(class = "section-label", "Stoichiometry Ranges"),
      uiOutput("stoich_controls"),
      
      hr(),
      
      # -- Generation mode --
      div(class = "section-label", "Generation Mode"),
      radioButtons("gen_mode", NULL,
                   choices = c(
                     "All combinations" = "all",
                     "Random subset" = "random",
                     "Quick homomer series" = "homomer",
                     "Custom ratios" = "custom"
                   ),
                   selected = "all"),
      conditionalPanel(
        "input.gen_mode == 'random'",
        numericInput("n_random", "How many random models?", 10, min = 1, max = 500)
      ),
      conditionalPanel(
        "input.gen_mode == 'custom'",
        textAreaInput("custom_ratios", "Enter ratios (one per line):",
                      placeholder = "2 1\n3 1\n1 1\n4 2",
                      rows = 5, width = "100%"),
        div(class = "help-text", "Space-separated copy numbers, one combo per line. Order matches your proteins.")
      ),
      
      br(),
      actionButton("generate_btn", "Generate Models", class = "btn-generate btn-block", width = "100%"),
      
      hr(),
      
      # -- Export format --
      div(class = "section-label", "Export Options"),
      checkboxGroupInput("export_fmt", NULL,
                         choices = c("AlphaFold3 JSON" = "json",
                                     "Multimer FASTA" = "fasta",
                                     "Summary Text" = "txt",
                                     "Manifest CSV" = "csv"),
                         selected = c("json", "fasta", "csv"))
    ),
    
    # ── Main Panel ──
    mainPanel(
      width = 9,
      
      uiOutput("stats_row"),
      
      tabsetPanel(
        id = "main_tabs", type = "tabs",
        
        # -- Tab: Sequences --
        tabPanel(
          "Sequences", value = "tab_seq",
          br(),
          DTOutput("protein_table"),
          uiOutput("sequence_details")
        ),
        
        # -- Tab: Models --
        tabPanel(
          "Models", value = "tab_models",
          br(),
          fluidRow(
            column(6,
                   actionButton("select_all_models", "Select All", class = "btn-sm btn-default"),
                   actionButton("deselect_all_models", "Deselect All", class = "btn-sm btn-default"),
                   actionButton("select_safe_models", "Select Safe (<3000 aa)", class = "btn-sm btn-info")
            ),
            column(6, style = "text-align: right;",
                   uiOutput("selection_count")
            )
          ),
          br(),
          DTOutput("combo_table")
        ),
        
        # -- Tab: Preview --
        tabPanel(
          "Preview", value = "tab_preview",
          br(),
          fluidRow(
            column(6, selectInput("preview_model", "Select model to preview:", choices = NULL, width = "100%")),
            column(6, br(), uiOutput("preview_model_info"))
          ),
          tabsetPanel(
            type = "pills",
            tabPanel("JSON", br(), verbatimTextOutput("json_preview")),
            tabPanel("FASTA", br(), verbatimTextOutput("fasta_preview")),
            tabPanel("Summary", br(), verbatimTextOutput("summary_preview"))
          )
        ),
        
        # -- Tab: Download --
        tabPanel(
          "Download", value = "tab_dl",
          br(),
          
          div(class = "dl-section",
              h4("Bulk Download (Selected Models)"),
              p("Downloads only the models you selected in the Models tab."),
              fluidRow(
                column(4, downloadButton("download_zip", "Download Selected (ZIP)", class = "btn-primary")),
                column(8, uiOutput("dl_summary"))
              )
          ),
          
          div(class = "dl-section",
              h4("Single Model Download"),
              fluidRow(
                column(5, selectInput("dl_single_model", "Pick a model:", choices = NULL, width = "100%")),
                column(7, br(),
                       downloadButton("download_single_json", "JSON", class = "btn-sm btn-default"),
                       downloadButton("download_single_fasta", "FASTA", class = "btn-sm btn-default"),
                       downloadButton("download_single_txt", "Summary", class = "btn-sm btn-default")
                )
              )
          ),
          
          div(class = "dl-section",
              h4("Full Manifest"),
              p("CSV of all generated models (selected or not) for your records."),
              downloadButton("download_manifest", "Download Manifest CSV", class = "btn-sm btn-default")
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    proteins = NULL,
    combos   = NULL,
    selected = c()
  )
  
  # ── Parse ──
  observeEvent(input$parse_btn, {
    text <- ""
    if (input$input_mode == "paste") {
      text <- input$fasta_text
    } else if (!is.null(input$fasta_file)) {
      text <- paste(readLines(input$fasta_file$datapath, warn = FALSE), collapse = "\n")
    }
    
    if (nchar(trimws(text)) == 0) {
      showNotification("Provide sequence data first.", type = "error")
      return()
    }
    
    proteins <- parse_fasta(text)
    if (is.null(proteins) || length(proteins) == 0) {
      showNotification("Could not parse any sequences. Check your FASTA format.", type = "error")
      return()
    }
    
    for (p in proteins) {
      err <- validate_sequence(p$sequence)
      if (!is.null(err)) {
        showNotification(paste0(p$id, ": ", err), type = "warning", duration = 8)
      }
    }
    
    rv$proteins <- proteins
    rv$combos <- NULL
    rv$selected <- c()
    showNotification(paste(length(proteins), "sequence(s) parsed successfully."), type = "message")
    updateTabsetPanel(session, "main_tabs", selected = "tab_seq")
  })
  
  output$parsed_summary_sidebar <- renderUI({
    req(rv$proteins)
    n <- length(rv$proteins)
    total_aa <- sum(sapply(rv$proteins, function(p) nchar(p$sequence)))
    tags$div(class = "protein-loaded", style = "margin-top: 10px;",
             tags$strong(paste0(n, " protein(s) loaded")),
             tags$br(),
             tags$span(style = "font-size: 0.85em; color: #555;",
                       paste0("Total: ", format(total_aa, big.mark = ","), " residues")),
             tags$br(),
             lapply(rv$proteins, function(p) {
               tags$span(class = "sequence-badge", paste0(p$id, " (", nchar(p$sequence), " aa)"))
             })
    )
  })
  
  # ── Protein Table ──
  output$protein_table <- renderDT({
    req(rv$proteins)
    df <- data.frame(
      ID = sapply(rv$proteins, function(p) p$id),
      `Length (aa)` = sapply(rv$proteins, function(p) nchar(p$sequence)),
      `First 40 residues` = sapply(rv$proteins, function(p) {
        s <- p$sequence
        if (nchar(s) > 40) paste0(substr(s, 1, 40), "...") else s
      }),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(pageLength = 20, dom = 'tip'), rownames = FALSE,
              selection = "single")
  })
  
  output$sequence_details <- renderUI({
    req(rv$proteins)
    sel <- input$protein_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      return(div(class = "help-text", style = "margin-top: 10px;",
                 "Click a row to view the full sequence."))
    }
    prot <- rv$proteins[[sel]]
    tags$div(style = "margin-top: 15px;",
             tags$h5(paste0(prot$id, " (", nchar(prot$sequence), " aa)")),
             tags$pre(style = "background: #1e1e1e; color: #7ec8e3; padding: 12px; border-radius: 8px;
                       font-family: Consolas, monospace; font-size: 0.85em; word-wrap: break-word; white-space: pre-wrap;",
                      prot$sequence)
    )
  })
  
  # ── Stoichiometry Controls ──
  output$stoich_controls <- renderUI({
    req(rv$proteins)
    controls <- lapply(seq_along(rv$proteins), function(i) {
      prot <- rv$proteins[[i]]
      div(class = "stoich-protein-block",
          tags$strong(prot$id),
          tags$span(style = "color: #888; font-size: 0.85em;", paste0("  (", nchar(prot$sequence), " aa)")),
          fluidRow(
            column(6, numericInput(paste0("min_", i), "Min", value = 1, min = 0, max = 10, step = 1)),
            column(6, numericInput(paste0("max_", i), "Max", value = 3, min = 1, max = 10, step = 1))
          )
      )
    })
    do.call(tagList, controls)
  })
  
  # ── Stats Row ──
  output$stats_row <- renderUI({
    n_prot <- if (!is.null(rv$proteins)) length(rv$proteins) else 0
    n_models <- if (!is.null(rv$combos)) length(rv$combos) else 0
    n_warn <- if (!is.null(rv$combos)) sum(sapply(rv$combos, function(c) c$warning)) else 0
    n_sel <- length(rv$selected)
    
    fluidRow(
      column(3, div(class = paste("stat-card", if (n_prot > 0) "success" else ""),
                    div(class = "stat-number", n_prot),
                    div(class = "stat-label", "Proteins Loaded"))),
      column(3, div(class = "stat-card",
                    div(class = "stat-number", n_models),
                    div(class = "stat-label", "Models Generated"))),
      column(3, div(class = paste("stat-card", if (n_warn > 0) "warning" else ""),
                    div(class = "stat-number", n_warn),
                    div(class = "stat-label", "Over 3000 aa"))),
      column(3, div(class = paste("stat-card", if (n_sel > 0) "success" else ""),
                    div(class = "stat-number", n_sel),
                    div(class = "stat-label", "Selected")))
    )
  })
  
  # ── Generate Combinations ──
  observeEvent(input$generate_btn, {
    req(rv$proteins)
    n_prot <- length(rv$proteins)
    prot_ids <- sapply(rv$proteins, function(p) p$id)
    
    ranges <- lapply(seq_len(n_prot), function(i) {
      mn <- input[[paste0("min_", i)]]
      mx <- input[[paste0("max_", i)]]
      if (is.null(mn)) mn <- 1
      if (is.null(mx)) mx <- 3
      mn <- max(0, mn)
      mx <- max(mn, mx)
      seq(mn, mx)
    })
    
    combos <- NULL
    
    if (input$gen_mode == "custom") {
      text <- input$custom_ratios
      if (is.null(text) || nchar(trimws(text)) == 0) {
        showNotification("Enter at least one ratio in the custom ratios box.", type = "error")
        return()
      }
      combos <- generate_custom_combinations(rv$proteins, text)
      
    } else if (input$gen_mode == "homomer") {
      combos <- list()
      for (i in seq_along(rv$proteins)) {
        for (n in ranges[[i]]) {
          counts <- rep(0, n_prot)
          counts[i] <- n
          names(counts) <- prot_ids
          total_res <- sum(counts * sapply(rv$proteins, function(p) nchar(p$sequence)))
          combos[[length(combos) + 1]] <- list(
            model_name  = format_model_name(counts, prot_ids),
            display     = format_ratio_display(counts, prot_ids),
            counts      = counts,
            total_chains = sum(counts),
            total_residues = total_res,
            warning     = total_res > 3000
          )
        }
      }
      
    } else if (input$gen_mode == "random") {
      n_want <- input$n_random
      if (is.null(n_want)) n_want <- 10
      combos <- generate_random_combinations(rv$proteins, ranges, n = n_want)
      
    } else {
      total_combos <- prod(sapply(ranges, length))
      if (total_combos > 500) {
        showNotification(
          paste0(format(total_combos, big.mark = ","), " combinations. Consider Random mode or narrower ranges."),
          type = "warning", duration = 8)
        if (total_combos > 5000) {
          showNotification("Blocked: >5,000 combinations. Use Random mode.", type = "error")
          return()
        }
      }
      combos <- generate_combinations(rv$proteins, ranges)
    }
    
    rv$combos <- combos
    rv$selected <- seq_along(combos)
    
    choices <- sapply(combos, function(c) c$model_name)
    display <- paste0(sapply(combos, function(c) c$display),
                      " (", sapply(combos, function(c) format(c$total_residues, big.mark = ",")), " aa)")
    names(choices) <- display
    updateSelectInput(session, "preview_model", choices = choices)
    updateSelectInput(session, "dl_single_model", choices = choices)
    
    n_warn <- sum(sapply(combos, function(c) c$warning))
    msg <- paste0(length(combos), " model(s) generated. All selected for download.")
    if (n_warn > 0) msg <- paste0(msg, " ", n_warn, " exceed 3,000 residues.")
    showNotification(msg, type = "message", duration = 6)
    
    updateTabsetPanel(session, "main_tabs", selected = "tab_models")
  })
  
  # ── Model Table with Selection ──
  output$combo_table <- renderDT({
    req(rv$combos)
    df <- data.frame(
      `#` = seq_along(rv$combos),
      Model = sapply(rv$combos, function(c) c$model_name),
      Stoichiometry = sapply(rv$combos, function(c) c$display),
      Chains = sapply(rv$combos, function(c) c$total_chains),
      `Total Residues` = sapply(rv$combos, function(c) format(c$total_residues, big.mark = ",")),
      Status = sapply(rv$combos, function(c) if (c$warning) "OVER LIMIT" else "OK"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    datatable(df,
              options = list(pageLength = 25, dom = 'ftip', scrollX = TRUE,
                             columnDefs = list(list(className = 'dt-center', targets = c(3, 4, 5)))),
              rownames = FALSE,
              selection = list(mode = "multiple", selected = rv$selected)) %>%
      formatStyle("Status",
                  color = styleEqual(c("OVER LIMIT", "OK"), c("#e74c3c", "#27ae60")),
                  fontWeight = "bold")
  })
  
  observeEvent(input$combo_table_rows_selected, {
    rv$selected <- input$combo_table_rows_selected
  }, ignoreNULL = FALSE)
  
  observeEvent(input$select_all_models, {
    req(rv$combos)
    proxy <- dataTableProxy("combo_table")
    selectRows(proxy, seq_along(rv$combos))
    rv$selected <- seq_along(rv$combos)
  })
  
  observeEvent(input$deselect_all_models, {
    proxy <- dataTableProxy("combo_table")
    selectRows(proxy, NULL)
    rv$selected <- c()
  })
  
  observeEvent(input$select_safe_models, {
    req(rv$combos)
    safe <- which(!sapply(rv$combos, function(c) c$warning))
    proxy <- dataTableProxy("combo_table")
    selectRows(proxy, safe)
    rv$selected <- safe
  })
  
  output$selection_count <- renderUI({
    n <- length(rv$selected)
    total <- if (!is.null(rv$combos)) length(rv$combos) else 0
    if (total > 0) {
      tags$span(style = "font-weight: 600; color: #0f3460;",
                paste0(n, " / ", total, " models selected"))
    }
  })
  
  # ── Preview ──
  selected_combo <- reactive({
    req(rv$combos, input$preview_model)
    idx <- which(sapply(rv$combos, function(c) c$model_name) == input$preview_model)
    if (length(idx) > 0) rv$combos[[idx[1]]] else NULL
  })
  
  output$preview_model_info <- renderUI({
    combo <- selected_combo()
    req(combo)
    tags$div(style = "font-size: 0.9em; padding: 8px; background: #f0f7ff; border-radius: 6px;",
             tags$strong("Stoichiometry: "), combo$display, tags$br(),
             tags$strong("Chains: "), combo$total_chains, " | ",
             tags$strong("Residues: "), format(combo$total_residues, big.mark = ","),
             if (combo$warning) tags$span(class = "warning-text", " -- Over 3,000 aa limit")
    )
  })
  
  output$json_preview <- renderText({
    req(selected_combo())
    combo_to_json(selected_combo(), rv$proteins)
  })
  
  output$fasta_preview <- renderText({
    req(selected_combo())
    combo_to_fasta(selected_combo(), rv$proteins)
  })
  
  output$summary_preview <- renderText({
    req(selected_combo())
    combo_to_summary(selected_combo(), rv$proteins)
  })
  
  # ── Downloads ──
  output$dl_summary <- renderUI({
    req(rv$combos)
    n_sel <- length(rv$selected)
    fmts <- input$export_fmt
    n_fmt <- length(fmts[fmts != "csv"])
    n_files <- n_sel * n_fmt
    tags$div(style = "margin-top: 8px; font-size: 0.9em;",
             paste0(n_sel, " models x ", n_fmt, " format(s) = ", n_files, " files + manifest"))
  })
  
  output$download_zip <- downloadHandler(
    filename = function() {
      paste0("AF_stoichiometry_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(rv$combos, rv$proteins)
      
      sel_indices <- rv$selected
      if (length(sel_indices) == 0) {
        showNotification("No models selected. Go to Models tab and select rows.", type = "error")
        # Write an empty zip so download handler doesn't crash
        tmpdir <- file.path(tempdir(), "empty_export")
        dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
        note_file <- file.path(tmpdir, "NO_MODELS_SELECTED.txt")
        writeLines("No models were selected for export.", note_file)
        zip::zip(file, files = "NO_MODELS_SELECTED.txt", root = tmpdir)
        unlink(tmpdir, recursive = TRUE)
        return()
      }
      
      selected_combos <- rv$combos[sel_indices]
      
      tmpdir <- file.path(tempdir(), paste0("af_export_", as.integer(Sys.time())))
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
      
      fmts <- input$export_fmt
      filenames <- c()
      manifest_rows <- list()
      
      for (combo in selected_combos) {
        name <- combo$model_name
        
        if ("json" %in% fmts) {
          fp <- file.path(tmpdir, paste0(name, ".json"))
          writeLines(combo_to_json(combo, rv$proteins), fp)
          filenames <- c(filenames, paste0(name, ".json"))
        }
        if ("fasta" %in% fmts) {
          fp <- file.path(tmpdir, paste0(name, ".fasta"))
          writeLines(combo_to_fasta(combo, rv$proteins), fp)
          filenames <- c(filenames, paste0(name, ".fasta"))
        }
        if ("txt" %in% fmts) {
          fp <- file.path(tmpdir, paste0(name, "_summary.txt"))
          writeLines(combo_to_summary(combo, rv$proteins), fp)
          filenames <- c(filenames, paste0(name, "_summary.txt"))
        }
        
        manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
          Model = name,
          Stoichiometry = combo$display,
          Total_Chains = combo$total_chains,
          Total_Residues = combo$total_residues,
          Warning = if (combo$warning) "OVER_3000" else "",
          stringsAsFactors = FALSE
        )
      }
      
      manifest <- do.call(rbind, manifest_rows)
      write.csv(manifest, file.path(tmpdir, "manifest.csv"), row.names = FALSE)
      filenames <- c(filenames, "manifest.csv")
      
      zip::zip(file, files = filenames, root = tmpdir)
      unlink(tmpdir, recursive = TRUE)
    },
    contentType = "application/zip"
  )
  
  dl_combo <- reactive({
    req(rv$combos, input$dl_single_model)
    idx <- which(sapply(rv$combos, function(c) c$model_name) == input$dl_single_model)
    if (length(idx) > 0) rv$combos[[idx[1]]] else NULL
  })
  
  output$download_single_json <- downloadHandler(
    filename = function() paste0(dl_combo()$model_name, ".json"),
    content = function(file) writeLines(combo_to_json(dl_combo(), rv$proteins), file)
  )
  
  output$download_single_fasta <- downloadHandler(
    filename = function() paste0(dl_combo()$model_name, ".fasta"),
    content = function(file) writeLines(combo_to_fasta(dl_combo(), rv$proteins), file)
  )
  
  output$download_single_txt <- downloadHandler(
    filename = function() paste0(dl_combo()$model_name, "_summary.txt"),
    content = function(file) writeLines(combo_to_summary(dl_combo(), rv$proteins), file)
  )
  
  output$download_manifest <- downloadHandler(
    filename = function() paste0("manifest_", format(Sys.time(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(rv$combos)
      rows <- lapply(rv$combos, function(combo) {
        data.frame(
          Model = combo$model_name,
          Stoichiometry = combo$display,
          Total_Chains = combo$total_chains,
          Total_Residues = combo$total_residues,
          Warning = if (combo$warning) "OVER_3000" else "",
          stringsAsFactors = FALSE
        )
      })
      write.csv(do.call(rbind, rows), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
