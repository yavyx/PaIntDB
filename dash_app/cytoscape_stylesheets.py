default = [
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
        'selector': ':selected',
        'style': {
            'border-width': 5,
            'border-color': 'black',
            'border-opacity': 1,
            'opacity': 1,
            'font-size': 12,
            'z-index': 9999
        }
    }
]

fold_change = default + \
    [
        {
            'selector': '[significanceSource = "TnSeq"]',
            'style': {
                'background-color': 'gray'
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
                'background-color': '#1c2fff'  # blue
            }
        },
    ]


combined = default + \
    [
        {
            'selector': '[significanceSource = "RNASeq"]',
            'style': {
                'background-color': '#26e81c'  # green
            }
        },
        {
            'selector': '[significanceSource = "TnSeq"]',
            'style': {
                'background-color': '#1c2fff'  # blue
            }
        },
        {
            'selector': '[significanceSource = "both"]',
            'style': {
                'background-color': 'red'
            }
        }
    ]

