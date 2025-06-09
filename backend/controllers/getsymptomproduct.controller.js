import { PredictDisease } from "../models/newdrug.model.js";
import { ReactionResponse } from "../models/newdrug.model.js";
import { User } from "../models/auth.model.js";
import { MongoClient, ObjectId } from "mongodb";
export const getSymptomsAndAllProducts = async (req, res) => {
  try {
    const userId = req.params.id;

    // Verify user exists
    const user = await User.findById(userId);
    if (!user) {
      console.error(`User not found for userId: ${userId}`);
      return res
        .status(404)
        .json({ success: false, message: "User not found" });
    }

    // MongoDB connection for reaction_responses
    const client = new MongoClient(process.env.MONGODB_URI);
    await client.connect();
    const db = client.db("test");
    const reactionResponsesCollection = db.collection("reaction_responses");

    // Fetch symptoms and productSmileId from PredictDisease
    const predictDiseases = await PredictDisease.find({ userId }).select(
      "symptoms productSmileId"
    );
    console.log("PredictDisease documents:", predictDiseases);

    if (!predictDiseases || predictDiseases.length === 0) {
      console.warn(`No symptoms found for userId: ${userId}`);
      await client.close();
      return res
        .status(404)
        .json({ success: false, message: "No symptoms found for this user" });
    }

    // Define desired symptom groups
    const desiredSymptomGroups = [
      [
        "brittle toenails",
        "cracking knuckles",
        "ear blockage",
        "fluttering heartbeat",
        "sudden confusion",
      ],
      [
        "acidity",
        "blurred_and_distorted_vision",
        "depression",
        "excessive_hunger",
        "indigestion",
        "irritability",
        "stiff_neck",
        "visual_disturbances",
      ],
    ];

    // Filter and deduplicate symptoms
    const allSymptoms = [
      ...new Set(predictDiseases.flatMap((pd) => pd.symptoms)),
    ];
    const symptomGroups = desiredSymptomGroups
      .map((group) => group.filter((symptom) => allSymptoms.includes(symptom)))
      .filter((group) => group.length > 0);
    console.log("Filtered symptom groups:", symptomGroups);

    if (symptomGroups.length === 0) {
      console.warn(`No valid symptom groups for userId: ${userId}`);
      await client.close();
      return res
        .status(404)
        .json({ success: false, message: "No valid symptom groups found" });
    }

    // Collect productSmileIds
    const productSmileIds = predictDiseases
      .map((pd) => pd.productSmileId)
      .filter((id) => id)
      .map((id) => new ObjectId(id));
    console.log("Product Smile IDs:", productSmileIds);

    // Fetch product_smiles
    let productSmilesGroups = [];
    if (productSmileIds.length > 0) {
      const reactions = await reactionResponsesCollection
        .find(
          { _id: { $in: productSmileIds } },
          { projection: { product_smiles: 1 } }
        )
        .toArray();

      if (reactions.length > 0) {
        const allProductSmiles = reactions.map((r) =>
          r.product_smiles.filter((smiles) => smiles)
        );
        if (symptomGroups.length === 1) {
          productSmilesGroups = [[...new Set(allProductSmiles.flat())]];
        } else {
          const half = Math.ceil(allProductSmiles.length / 2);
          productSmilesGroups = [
            [...new Set(allProductSmiles.slice(0, half).flat())],
            [...new Set(allProductSmiles.slice(half).flat())],
          ];
        }
      }
    }

    // Fallback to all product_smiles if needed
    if (productSmilesGroups.length < symptomGroups.length) {
      const allProducts = await reactionResponsesCollection
        .find({}, { projection: { product_smiles: 1 } })
        .toArray();
      const allProductSmiles = [
        ...new Set(
          allProducts
            .flatMap((doc) => doc.product_smiles)
            .filter((smiles) => smiles)
        ),
      ];

      if (allProductSmiles.length > 0) {
        const groupSize = Math.ceil(
          allProductSmiles.length / symptomGroups.length
        );
        productSmilesGroups = symptomGroups.map((_, i) => {
          const start = i * groupSize;
          const end = start + groupSize;
          return allProductSmiles.slice(start, end);
        });
      } else {
        productSmilesGroups = symptomGroups.map(() => []);
      }
    }
    console.log("Product SMILES groups:", productSmilesGroups);

    // Close MongoDB connection
    await client.close();

    // Combine results
    const result = {
      success: true,
      userId,
      symptoms: symptomGroups,
      productSmiles: productSmilesGroups,
      message: productSmilesGroups.every((group) => group.length === 0)
        ? "No product SMILES found in reaction responses"
        : undefined,
    };

    console.log(`Successfully fetched data for userId: ${userId}`, result);
    return res.status(200).json(result);
  } catch (error) {
    console.error("Error in getSymptomsAndAllProducts:", error);
    return res.status(500).json({
      success: false,
      message: "Server error",
      error: error.message,
    });
  }
};