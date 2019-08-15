// NOTE: See this: https://github.com/brentp/slivar#how-it-works

// NOTE: Sets up default values. Since ducktape is using basic JS, we can't use
// ES5 defaults in functions

function allelic_balance_high_quality(sample) {
  // This function ensures that the allelic balance observed
  // is what is expected for given zygosity
  if (sample.alts == 0) return sample.AB <= 0.01
  else if (sample.alts == 1) return sample.AB >= 0.2 && sample.AB <= 0.8
  return sample.AB >= 0.9
}

function high_quality(sample, depth, gq) {
  if (depth == undefined) depth = 5
  if (gq == undefined) gq = 30
  return sample.DP > depth && sample.GQ > gq
}

function uniparental_disomy(proband, mom, dad) {
  if (mom == undefined) return false
  if (dad == undefined) return false
  return (
    (proband.alts == 0 || proband.alts == 2) && ((mom.alts == 2 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 2))
  )
}

function sample_meets_dominant(sample) {
  if (sample.affected == true && sample.alts == 1) return true
  if (sample.affected == false && sample.alts == 0) return true
  return false
}

function dominant(proband, mom, dad) {
  if (mom == undefined) mom = {}
  if (dad == undefined) dad = {}

  if (mom.affected == false && dad.affected == false) return false

  if (proband.alts == 1) {
    if (sample_meets_dominant(mom) == false) return false
    if (sample_meets_dominant(dad) == false) return false
  } else {
    return false
  }

  return true
}

function denovo(proband, mom, dad) {
  if (mom == undefined) return false
  if (dad == undefined) return false

  if (mom.affected || dad.affected) return false

  if (proband.alts == 1) {
    if (mom.alts != undefined && mom.alts != 0) return false
    if (dad.alts != undefined && dad.alts != 0) return false
  } else {
    return false
  }

  return true
}

function x_linked_denovo(proband, mom, dad) {
  if (mom == undefined) mom = {}
  if (dad == undefined) dad = {}

  if (mom.affected || dad.affected) return false

  if (proband.alts >= 1) {
    if (proband.sex != undefined && proband.sex != 'male') return false
    if (mom.alts != undefined && mom.alts != 0) return false
    if (dad.alts != undefined && dad.alts != 0) return false
  } else {
    return false
  }

  return true
}

function homozygous_recessive(proband, mom, dad) {
  if (mom == undefined) mom = {}
  if (dad == undefined) dad = {}

  if (mom.affected || dad.affected) return false

  if (proband.alts == 2) {
    if (mom.alts != undefined && mom.alts != 1) return false
    if (dad.alts != undefined && dad.alts != 1) return false
  } else {
    return false
  }

  return true
}

function x_linked_homozygous_recessive(proband, mom, dad) {
  if (mom == undefined) mom = {}
  if (dad == undefined) dad = {}

  if (mom.affected || dad.affected) return false

  if (proband.alts >= 1) {
    if (proband.sex != undefined && proband.sex != 'male') return false
    if (mom.alts != undefined && mom.alts == 0) return false
    if (dad.alts != undefined && dad.alts != 0) return false
  } else {
    return false
  }

  return true
}

function compound_heterozygous_side(proband, mom, dad) {
  // NOTE: this identifies if a variant is potentially a compound heterozygous
  // NOTE: This is the 1st half, it then needs to run through slivar comphets to do in-phase transmission

  if (mom == undefined) mom = {}
  if (dad == undefined) dad = {}

  if (mom.affected || dad.affected) return false

  if (proband.alts == 1) {
    if (mom.alts != undefined && mom.alts == 2) return false
    if (dad.alts != undefined && dad.alts == 2) return false
  } else {
    return false
  }

  return true
}

function remove_compound_heterozygous_side(INFO, key) {
  // NOTE: this is useful for removing all the intermediate comphet side variants
  // after running slivar comp-hets

  if (key == undefined) key = 'compound_heterozygous_side'

  if (key in INFO) return false
  else return true
}

function proband_has_variant(proband) {
  return proband.alts > 0
}

function present_in_database(INFO, key, position) {
  if (key == undefined) key = 'CSQ'
  if (position == undefined) position = 1
  if (!(key in INFO)) return false

  var result = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')
    if (data[position]) result = true
    else result = false
  })

  return result
}

function includes_filter(INFO, key, position, types) {
  if (key == undefined) key = 'CSQ'
  if (position == undefined) position = 1
  if (types == undefined) types = []
  if (!(key in INFO)) return false

  var result = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')
    const element = data[position]
    types.forEach(function(type) {
      if (element.includes(type)) result = true
    })
  })

  return result
}

function match_filter(INFO, key, position, types) {
  if (key == undefined) key = 'CSQ'
  if (position == undefined) position = 1
  if (types == undefined) types = []
  if (!(key in INFO)) return false

  var result = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')
    const element = data[position]
    types.forEach(function(type) {
      if (element == type) result = true
    })
  })

  return result
}

function gte_filter(INFO, key, position, cutoff) {
  if (key == undefined) key = 'CSQ'
  if (position == undefined) position = 1
  if (cutoff == undefined) cutoff = 0.01
  if (!(key in INFO)) return false

  var impact = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')
    if (!data[position]) impact = false
    else {
      const scores = []
      data[position].split('&').forEach(function(score) {
        scores.push(parseFloat(score))
      })
      const score = Math.max.apply(null, scores)
      if (score >= parseFloat(cutoff)) impact = true
    }
  })

  return impact
}

function lte_filter(INFO, key, position, cutoff) {
  if (key == undefined) key = 'CSQ'
  if (position == undefined) position = 1
  if (cutoff == undefined) cutoff = 0.01
  if (!(key in INFO)) return false

  var impact = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')
    if (!data[position]) impact = false
    else {
      const scores = []
      data[position].split('&').forEach(function(score) {
        scores.push(parseFloat(score))
      })
      const score = Math.min.apply(null, scores)
      if (score <= parseFloat(cutoff)) impact = true
    }
  })

  return impact
}

function maf_filter(INFO, key, position, cutoff) {
  if (key == undefined) key = 'CSQ'
  if (position == undefined) position = 1
  if (cutoff == undefined) cutoff = 0.01
  if (!(key in INFO)) return false

  var impact = false
  const effects = INFO[key].split(',')
  effects.forEach(function(effect) {
    const data = effect.split('|')
    if (!data[position]) impact = true
    else if (parseFloat(data[position]) <= parseFloat(cutoff)) impact = true
  })

  return impact
}

module.exports = {
  allelic_balance_high_quality,
  high_quality,
  uniparental_disomy,
  sample_meets_dominant,
  dominant,
  denovo,
  x_linked_denovo,
  homozygous_recessive,
  x_linked_homozygous_recessive,
  compound_heterozygous_side,
  remove_compound_heterozygous_side,
  proband_has_variant,
  present_in_database,
  includes_filter,
  match_filter,
  gte_filter,
  lte_filter,
  maf_filter,
}
