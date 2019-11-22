combined = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)',
            'border-width': 1,
            'border-color': 'black',
            'min-zoomed-font-size': 3,
            'min-zoomed-font-size': 10
        }
    },
    {
        'selector': '[significanceSource = "RNASeq"]',
        'style': {
            'background-color': 'green'
        }
    },
    {
        'selector': '[significanceSource = "TnSeq"]',
        'style': {
            'background-color': 'blue'
        }
    },
    {
        'selector': '[significanceSource = "both"]',
        'style': {
            'background-color': 'red'
        }
    }
]

fold_change = [
    {
        'selector': 'node',
        'style': {
            'content': 'data(label)',
            'border-width': 1,
            'border-color': 'black',
            'min-zoomed-font-size': 10
        }
    },
    {
        'selector': '[significanceSource = "TnSeq"]',
        'style': {
            'background-color': 'green'
        }
    },
    {
        'selector': '[log2FoldChange < 0]',
        'style': {
            'background-color': 'red'
        }
    },
    {
        'selector': '[log2FoldChange > 0]',
        'style': {
            'background-color': 'blue'
        }
    },
]
