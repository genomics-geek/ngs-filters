// NOTE: See this: https://github.com/brentp/slivar#how-it-works

function allelic_balance_high_quality(sample) {
  // This function ensures that the allelic balance observed
  // is what is expected for given zygosity
  if (sample.alts == 0) return sample.AB <= 0.01
  if (sample.alts == 1) return sample.AB >= 0.2 && sample.AB <= 0.8
  if (sample.alts == 2) return sample.AB >= 0.99
  return false
}

function high_quality(sample, depth, gq) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (depth == undefined) {
    depth = 5
  }
  if (gq == undefined) {
    gq = 30
  }

  return sample.DP > depth && sample.GQ > gq
}

function uniparental_disomy(proband, mom, dad) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (mom == undefined) {
    mom = {}
  }
  if (dad == undefined) {
    dad = {}
  }

  return (
    (proband.alts == 0 || proband.alts == 2) && ((mom.alts == 2 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 2))
  )
}

function denovo(proband, mom, dad) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (mom == undefined) {
    mom = {}
  }
  if (dad == undefined) {
    dad = {}
  }

  if (proband.alts == 1) {
    if (mom.alts != undefined && mom.alts != 0) {
      return false
    }
    if (dad.alts != undefined && dad.alts != 0) {
      return false
    }
  } else {
    return false
  }

  return true
}

function x_linked_denovo(proband, mom, dad) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (mom == undefined) {
    mom = {}
  }
  if (dad == undefined) {
    dad = {}
  }

  if (proband.alts >= 1) {
    if (proband.sex != undefined && proband.sex != 'male') {
      return false
    }
    if (mom.alts != undefined && mom.alts != 0) {
      return false
    }
    if (dad.alts != undefined && dad.alts != 0) {
      return false
    }
  } else {
    return false
  }

  return true
}

function homozygous_recessive(proband, mom, dad) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (mom == undefined) {
    mom = {}
  }
  if (dad == undefined) {
    dad = {}
  }

  if (proband.alts == 2) {
    if (mom.alts != undefined && mom.alts == 0) {
      return false
    }
    if (dad.alts != undefined && dad.alts == 0) {
      return false
    }
  } else {
    return false
  }

  return true
}

function x_linked_homozygous_recessive(proband, mom, dad) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (mom == undefined) {
    mom = {}
  }
  if (dad == undefined) {
    dad = {}
  }

  if (proband.alts >= 1) {
    if (proband.sex != undefined && proband.sex != 'male') {
      return false
    }
    if (mom.alts != undefined && mom.alts == 0) {
      return false
    }
    if (dad.alts != undefined && dad.alts != 0) {
      return false
    }
  } else {
    return false
  }

  return true
}

function compound_heterozygous_side(proband, mom, dad) {
  // NOTE: this identifies if a variant is potentially a compound heterozygous
  // NOTE: This is the 1st half, it then needs to run through slivar comphets to do in-phase transmission

  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (mom == undefined) {
    mom = {}
  }
  if (dad == undefined) {
    dad = {}
  }

  if (proband.alts == 1) {
    if (mom.alts != undefined && mom.alts == 0) {
      return false
    }
    if (dad.alts != undefined && dad.alts == 0) {
      return false
    }
  } else {
    return false
  }

  return true
}

function remove_compound_heterozygous_side(INFO, key) {
  // NOTE: this is useful for removing all the intermediate comphet side variants
  // after running slivar comp-hets

  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (key == undefined) key = 'compound_heterozygous_side'

  if (key in INFO) {
    return false
  } else {
    return true
  }
}

function in_hgmd(INFO, key, position) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (key == undefined) {
    key = 'CSQ'
  }
  if (position == undefined) {
    position = 28
  }
  if (!(key in INFO)) {
    return false
  }

  var impact = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')
    if (data[position]) {
      impact = true
    }
  })

  return impact
}

function nonsynonymous(INFO, key, position) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (key == undefined) {
    key = 'CSQ'
  }
  if (position == undefined) {
    position = 1
  }
  if (!(key in INFO)) {
    return false
  }

  const types = ['missense', 'splice', 'insertion', 'deletion', 'frameshift', 'stop', 'start', 'coding sequence']

  var impact = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')
    const consequence = data[position]
    types.forEach(function(type) {
      if (consequence.includes(type)) impact = true
    })
  })

  return impact
}

function cohort_af(INFO, key, position, cutoff) {
  // NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
  // ES5 defaults in functions
  if (key == undefined) {
    key = 'CSQ'
  }
  if (position == undefined) {
    position = 35
  }
  if (cutoff == undefined) {
    cutoff = 0.01
  }
  if (!(key in INFO)) {
    return false
  }

  var impact = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')

    if (data[position] == undefined) {
      imppact = true
    } else if (data[position] <= cutoff) {
      impact = true
    }
  })

  return impact
}

module.exports = {
  allelic_balance_high_quality,
  high_quality,
  uniparental_disomy,
  denovo,
  x_linked_denovo,
  homozygous_recessive,
  x_linked_homozygous_recessive,
  compound_heterozygous_side,
  remove_compound_heterozygous_side,
  in_hgmd,
  nonsynonymous,
  cohort_af,
}
